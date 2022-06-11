//
//  main.cpp
//  essentiaCPP_test1
//
//  Created by Nicholas Solem on 10/28/21.
//  Copyright © 2021 Nicholas Solem. All rights reserved.
//
#include<stdio.h>
#include<stdlib.h>
#include <cfloat> // for FLT_EPSILON
#include <iostream>
#include <dirent.h>
//#include <filesystem>
#include "essentia.h"
#include "pool.h"
#include "algorithmfactory.h"
#include "essentiamath.h" // includes <algorithm>, <vector>, <numeric> (needed by sort_indexes)
#include "noteToHertzMap.h"

#define MEASURE_TIME 1
#define UIOWA 1
#if MEASURE_TIME
#include <ctime>
#endif

// make array of filenames
// and it could spit out the indices to filenames as well as their respective f0

std::vector<std::string> parseUIOWAFile(std::string fnToParse){
    size_t pos = 0;
    std::vector<std::string> token;
    token.reserve(4);
    std::string delimiter = ".";
    int count = 0;
    std::cout << fnToParse << "\n";
    while ((pos = fnToParse.find(delimiter)) != std::string::npos) {
//        token[count] = fnToParse.substr(0, pos);
        token.push_back(fnToParse.substr(0, pos));
        //std::cout << token[count] << std::endl;
        fnToParse.erase(0, pos + delimiter.length());
        count++;
    }
    return token;
}
// function to return indices of sorted vector without changing that vector
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {
// this function was found on stackoverflow
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}
void a2db(essentia::Real& amplitude){
    amplitude = 20.f * log10f(amplitude);
}

void d_print(const char path[], int N, essentia::Real* array){
    // "/Users/nicholassolem/development/essentiaTestFile.txt"
    FILE *fp = fopen(path, "w");
    if (fp != NULL) {
        for (int i = 0; i < N; i++) fprintf(fp, "%f\n", array[i]);
        fclose(fp);
    }
    else
        printf("file NULL");
}

using namespace essentia::standard;

int main(int argc, char* argv[]) {
    makeMap(); // build noteToHertzMap
    const int sampleRate = 44100;
    const int frameSize = 512;
    const int frameSizeBig = frameSize * 4;
    const int hopSize = 2000; // WAY TOO BIG OF ILE STILL. perhaps will be necessary to make adaptive hop size, which is shorter near transients?
    const float wave_f0 = (float)sampleRate / (float)frameSize;
    const float lowest_f0 = (float)sampleRate / (float)frameSizeBig;
    const float specNormalizerBig = 1.f / (float)frameSize;
    const int N_HARMONICS = 64;
    const int N_OUT_FFT = 128; // next power of 2
    
    std::vector<essentia::Real> startTimes, endTimes;
    
    essentia::init();
    AlgorithmFactory& factory = essentia::standard::AlgorithmFactory::instance();

    std::string initFn = "/Users/nicholassolem/development/Yin/saw100Hz.wav";// "/Users/nicholassolem/Local/master samples/violin_loops/Royer1.wav";
//    initFn = "/Users/nicholassolem/Documents/MATLAB/wtsource/uiowa/piano/Piano.pp.Db1.wav";
    Algorithm* audio    = factory.create("MonoLoader", "filename", initFn, "sampleRate", sampleRate);

    Algorithm* fc       = factory.create("FrameCutter",
                                         "frameSize", frameSize,
                                         "hopSize", hopSize);
    Algorithm* fc_big   = factory.create("FrameCutter",
                                         "frameSize", frameSizeBig,
                                         "hopSize", hopSize); // need to confirm centering of windows
//    Algorithm* slicer= factory.create("Slicer", "endTimes", endTimes, "sampleRate", sampleRate, "startTimes", startTimes, "timeUnits", "samples");
    Algorithm* w        = factory.create("Windowing",
                                      "type", "blackmanharris92");
    Algorithm* w_big    = factory.create("Windowing",
                                      "type", "blackmanharris92");
    Algorithm* dcBlk    = factory.create("DCRemoval", "sampleRate", sampleRate, "cutoffFrequency", lowest_f0); //without, there can be DC peaks which throws error for harmonic peaks
//    Algorithm* dcBlkBig = factory.create("DCRemoval", "sampleRate", sampleRate, "cutoffFrequency", 20.f); //without, there can be DC peaks which throws error for harmonic peaks
    Algorithm* spec     = factory.create("Spectrum");
    Algorithm* powSpec  = factory.create("PowerSpectrum");
//    Algorithm* mfcc  = factory.create("MFCC","sampleRate", sampleRate, "highFrequencyBound", sampleRate / 2); // could also just use this but throw out all but the first 20-30 of them
//    Algorithm* crepe  = factory.create("PitchCREPE", "batchSize", 64);
    Algorithm* bfcc     = factory.create("BFCC","sampleRate", sampleRate, "type", "power", "logType", "dbpow", "highFrequencyBound", sampleRate / 2, "inputSize", frameSizeBig / 2 + 1, "numberBands", 23, "numberCoefficients", 23);// up to 7700 Hz; if using fs>=44100, use 27 or 28 bands
//    Algorithm* barks    = factory.create("BarkBands", "sampleRate", sampleRate, "numberBands", 23); // up to 7700 Hz; if using fs>=44100, use 27 or 28 bands
    Algorithm* rms      = factory.create("RMS");
    Algorithm* loud     = factory.create("Loudness");
    Algorithm* pyp      = factory.create("PitchYinProbabilistic", "frameSize", frameSize, "hopSize", hopSize, "lowRMSThreshold", 0.1, "outputUnvoiced", "abs", "preciseTime", false, "sampleRate", sampleRate);
    Algorithm* sCmplxty = factory.create("SpectralComplexity", "magnitudeThreshold", 0.0055, "sampleRate", sampleRate);
    Algorithm* sPeaks   = factory.create("SpectralPeaks","sampleRate", sampleRate, "minFrequency", lowest_f0);
    Algorithm* hPeaks   = factory.create("HarmonicPeaks");
    Algorithm* inharm   = factory.create("Inharmonicity");
    Algorithm* flat     = factory.create("Flatness");
    
    Algorithm* rs       = factory.create("Resample", "inputSampleRate", sampleRate, "outputSampleRate", sampleRate, "quality", 0); // inputSampleRate will be overwritten each iteration. Quality=0 is best quality.
    Algorithm* rs_win   = factory.create("Windowing", "type", "hann");
    Algorithm* rs_spec  = factory.create("Spectrum");
    
    Algorithm* una_abs  = factory.create("UnaryOperatorStream", "type", "abs");
    Algorithm* min      = factory.create("MinMax", "type", "min");
    
    /////////// CONNECTING THE ALGORITHMS ////////////////
    std::cout << "-------- connecting algos ---------" << std::endl;
    // register the algorithms in the factory(ies)
    essentia::Pool pool;
    
    // Audio -> FrameCutter
    std::vector<essentia::Real> audioBuffer;

    audio->output("audio").set(audioBuffer);
//    slicer->input("audio").set(audioBuffer);   // frame creator

    // FrameCutter -> Windowing -> Spectrum
    std::vector<std::vector<essentia::Real>> framesMany;//, crepeActivations;
    std::vector<essentia::Real> frameSingle, frameSingleBig,frameSingleBigAbs, windowedFrame,windowedFrameBig,hpWinFrame, rsFrame, rsSpec, rsWin;
    std::vector<essentia::Real> spectrum, powSpectrum, barkCoeffs, bfccCoeffs, bfccBands,sFreqs,sMags,hFreqs,hMags;
    std::vector<essentia::Real> pitchPYP, voicedProbPYP;
        //,crepeTimes, crepeFreqs, crepeConf;
    essentia::Real loudnesses, rmss, pitches, complexities, inharmonicities, flatnesses, frameMinVal;
    int frameMinIdx;
    
    fc->input("signal").set(audioBuffer);   // frame creator
    fc->output("frame").set(frameSingle);
    fc_big->input("signal").set(audioBuffer);
    fc_big->output("frame").set(frameSingleBig);
    pyp->input("signal").set(audioBuffer);
    pyp->output("pitch").set(pitchPYP);
    pyp->output("voicedProbabilities").set(voicedProbPYP);
//    slicer->output("frame").set(framesMany);
    w->input("frame").set(frameSingle);
    w->output("frame").set(windowedFrame);
    w_big->input("frame").set(frameSingleBig);
    w_big->output("frame").set(windowedFrameBig);
    //DC remove
    dcBlk->input("signal").set(windowedFrameBig);
    dcBlk->output("signal").set(hpWinFrame);
    
    spec->input("frame").set(hpWinFrame);
    spec->output("spectrum").set(spectrum);
    
    powSpec->input("signal").set(hpWinFrame);
    powSpec->output("powerSpectrum").set(powSpectrum);

//    crepe->input("audio").set(audioBuffer);
//    crepe->output("time").set(crepeTimes);
//    crepe->output("frequency").set(crepeFreqs);
//    crepe->output("confidence").set(crepeConf);
//    crepe->output("activations").set(crepeActivations);
//
    rms->input("array").set(hpWinFrame /*frameSingle*/);
    rms->output("rms").set(rmss);
    loud->input("signal").set(hpWinFrame);
    loud->output("loudness").set(loudnesses);
    
//    barks->input("spectrum").set(powSpectrum);
//    barks->output("bands").set(barkCoeffs);

    bfcc->input("spectrum").set(powSpectrum);
    bfcc->output("bands").set(bfccBands);
    bfcc->output("bfcc").set(bfccCoeffs);

    sCmplxty->input("spectrum").set(spectrum);
    sCmplxty->output("spectralComplexity").set(complexities);
    
    sPeaks->input("spectrum").set(spectrum);    // "recommended that the input "spectrum" be computed by the Spectrum algorithm"
    sPeaks->output("frequencies").set(sFreqs);
    sPeaks->output("magnitudes").set(sMags);

//    pitches = f0;
    hPeaks->input("frequencies").set(sFreqs);   //need to set sPeaks out first
    hPeaks->input("magnitudes").set(sMags);
    hPeaks->input("pitch").set(pitches);        // takes in a REFERENCE so pitch can change with each iteration
    hPeaks->output("harmonicFrequencies").set(hFreqs);
    hPeaks->output("harmonicMagnitudes").set(hMags);
    
    inharm->input("frequencies").set(hFreqs);
    inharm->input("magnitudes").set(hMags);
    inharm->output("inharmonicity").set(inharmonicities);
    
    flat->input("array").set(spectrum);
    flat->output("flatness").set(flatnesses);
    
    // to get zero cross as 0th sample; might not matter consistently since we are discarding phase
    una_abs->input("array").set(frameSingleBig);
    una_abs->output("array").set(frameSingleBigAbs);
    min->input("array").set(frameSingleBigAbs);
    min->output("real").set(frameMinVal);
    min->output("int").set(frameMinIdx);
    // resample from big window to next pow 2 of bin truncation (so 128)
    rs->input("signal").set(frameSingleBig);
    rs->output("signal").set(rsFrame);
    rs_win->input("frame").set(rsFrame);
    rs_win->output("frame").set(rsWin);
    rs_spec->input("frame").set(rsWin);
    rs_spec->output("spectrum").set(rsSpec);
    
//    std::string path = "/Users/nicholassolem/development/audio for analysis"; //initFn;
//    std::string audioFilename =  "/Users/nicholassolem/development/audio for analysis/Piano.ff.Db7.wav";
//    std::string audioFilename2 = "/Users/nicholassolem/development/audio for analysis/Piano.mf.D4.wav";
    std::vector<std::string> audioFiles;// = {audioFilename, audioFilename2};
    const char *path = "/Users/nicholassolem/development/audio for analysis/plucked";
//    audioFiles.push_back("/Users/nicholassolem/development/audio for analysis/smooth.wav");
    DIR *d;
        struct dirent *dir;
        d = opendir(path);
        if (d) {
          while ((dir = readdir(d)) != NULL) {
              char destPath[500];
              char destName[120];
              strcpy(destPath, path);
              strcpy(destName, dir->d_name);
              strcat(destPath, "/");
              strcat(destPath, destName);
    //          printf("%s\n", destPath);
              int endOfName = (int)strlen(destPath);
              char type[4] = {destPath[endOfName-3], destPath[endOfName-2], destPath[endOfName-1] };
              if (strcmp(type, "wav") == 0){ // no difference, so it ends with 'wav'
                  audioFiles.push_back(destPath);
              }
          }
          closedir(d);
        for(std::string s : audioFiles)
            std::cout << s << "\n";
        }
    std::string outputFilename(path);
    outputFilename.append("/pluck_essentia.yaml");
    
    /////////// STARTING THE ALGORITHMS //////////////////
    std::cout << "-------- start processing " << path << " --------" << std::endl;
    int numFiles = (int)audioFiles.size();

    for (int n=0; n< numFiles; n++){
#if UIOWA
        std::vector<std::string> token = parseUIOWAFile(audioFiles[n]);
        
        std::string dynamic = token[1];
//        std::string note = token[2];
//        float f0 = noteToHertzMap[token[2]];
        float f0 = noteToHertzMap[token[4]];
        std::cout << f0 << std::endl;
        pitches = f0; // / (essentia::Real)sampleRate;
#endif
        audio = factory.create("MonoLoader", "filename", audioFiles[n], "sampleRate", sampleRate);
        audio->output("audio").set(audioBuffer);
        
        // create start & end times for slicer
        // use audioBuffer?
        int m;
//        size_t length = audioBuffer.size();
//        size_t startsRange = length - frameSize;
//        for (m = 0; m < 1; m++){
//            size_t start = (m * startsRange) / 1;
//            startTimes.push_back((essentia::Real)start);
//            size_t end = start + frameSize - 1;
//            endTimes.push_back((essentia::Real)end);
//        }
        audio->compute();
        fc       = factory.create("FrameCutter",
                                             "frameSize", frameSize,
                                             "hopSize", hopSize);
        fc_big   = factory.create("FrameCutter",
                                             "frameSize", frameSizeBig,
                                             "hopSize", hopSize);
        fc->input("signal").set(audioBuffer);   // frame creator
        fc->output("frame").set(frameSingle);
        fc_big->input("signal").set(audioBuffer);
        fc_big->output("frame").set(frameSingleBig);
#if MEASURE_TIME
        std::clock_t c_start = std::clock();
#endif
        pyp->compute();
#if MEASURE_TIME
        std::clock_t c_end = std::clock();
        long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
        std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
#endif
        m = 0;
        while (true) {
#if !UIOWA
            float f0 = pitchPYP[m];
            pitches = f0;
#else
            
#endif
            // compute a frame
            fc_big->compute();
            fc->compute();
//            frameSingle = framesMany[m++];
//            w->input("frame").set(frameSingle);
//            w_big->input("frame").set(frameSingleBig);
            // if it was the last one (ie: it was empty), then we're done.
            if (!(frameSingle.size()) || !(frameSingleBig.size())) {break;}

            // if the frame is silent, just drop it and go on processing
            bool smallSilent = essentia::isSilent(frameSingle);
            bool bigSilent = essentia::isSilent(frameSingleBig);
            if (smallSilent || bigSilent) continue;

            w->compute();
            w_big->compute();
            dcBlk->compute();
            spec->compute();
            powSpec->compute();
//            crepe->compute();
//            mfcc->compute();
            bfcc->compute();
//            barks->compute();
            rms->compute();
            loud->compute();
            sCmplxty->compute();
            sPeaks->compute();
            // convert magnitude to dB
            
            std::for_each(std::begin(sMags), std::end(sMags), [specNormalizerBig](essentia::Real& amplitude) {amplitude = 20.f * log10f((amplitude + FLT_EPSILON) * frameSize);});

            hPeaks->compute();
            inharm->compute();
            flat->compute();
            
            // get center index
            int c_idx = frameSizeBig / 2;
            // get wavelength
            essentia::Real wavelength = (essentia::Real)sampleRate / f0;
            
            una_abs->compute();
            min->compute();         // computes frameMinIdx
            int st_idx = c_idx - wavelength;
            int end_idx = c_idx + wavelength;
            // problem with doing this way: prioritizes zeros at beginning of frame, which makes the 0th frame all 0's due to zero-padding
//            int st_idx = frameMinIdx <= c_idx ? frameMinIdx : c_idx; // sample closest to 0 as long as it's before center of frame
//            int end_idx = st_idx + 2 * wavelength;
            std::vector<essentia::Real>::const_iterator first = frameSingleBig.begin() + st_idx;
            std::vector<essentia::Real>::const_iterator last = frameSingleBig.begin() + end_idx;
            std::vector<essentia::Real> wave(first, last);

            rs = factory.create("Resample", "inputSampleRate", (essentia::Real)wave.size(), "outputSampleRate", 512, "quality", 0);
            rs->input("signal").set(wave);
            rs->output("signal").set(rsFrame);
            rs->compute();
//            d_print("/Users/nicholassolem/development/essentia_wave.txt", (int)wave.size(), &wave[0]);
//            d_print("/Users/nicholassolem/development/essentia_rsFrame.txt", (int)rsFrame.size(), &rsFrame[0]);
            if (rsFrame.size() % 2 != 0)
                rsFrame.push_back(0.0);
            rs_win->compute();
            rs_spec->compute();
            rsSpec.resize((size_t)N_HARMONICS);
            for (int i = 0; i <rsSpec.size(); i++){
                rsSpec[i] = powf(rsSpec[i], 0.3f);
            }
//            d_print("/Users/nicholassolem/development/essentia_rsSpec.txt", frameSizeBig, &rsSpec[0]);
            pool.add("lowlevel.freqs", pitches / (essentia::Real)sampleRate);
//            pool.add("lowlevel.mfcc", mfccCoeffs);
            pool.add("lowlevel.bfcc", bfccCoeffs);
            pool.add("lowlevel.bark", bfccBands);
            pool.add("lowlevel.rms", rmss);
            pool.add("lowlevel.loudness", loudnesses);
//            pool.add("lowlevel.spectralFrequencies",sFreqs);
//            pool.add("lowlevel.spectralMagnitudes",sMags);
//            pool.add("lowlevel.harmonicFrequencies",hFreqs);
//            pool.add("lowlevel.harmonicMagnitudes", hMags);
            pool.add("lowlevel.spectralComplexity", complexities);
            pool.add("lowlevel.inharmonicity",inharmonicities);
//            pool.add("lowlevel.flatness",flatnesses);
            pool.add("lowlevel.resampledSpectrum", rsSpec);
//            pool.add("lowlevel.zeroCrossingRate",zcrs);
            m++;
        }
    }

    
    // CREATE SORTED REFERENCE TO INDEPENDENT FEATURE VECTORS TO ALLOW FOR CULLING iTH ELEMENT OF ALL
    // or should i just delete from pool?
    
    essentia::Pool cullPool;
    std::vector<essentia::Real> loudPooled, loudCull, rmsPooled, rmsCull, pitchPooled, pitchCull, cplxtyPooled, cplxtyCull, inharmPooled, inharmCull;
    std::vector<std::vector<essentia::Real>> barkPooled, barkCull, bfccPooled, bfccCull, rsSpecPooled, rsSpecCull;
    essentia::Real percentileLow = 0.25;
    essentia::Real percentileHigh = 0.99;
    loudPooled   = pool.value<std::vector<essentia::Real>>("lowlevel.loudness");
    rmsPooled    = pool.value<std::vector<essentia::Real>>("lowlevel.rms");
    pitchPooled  = pool.value<std::vector<essentia::Real>>("lowlevel.freqs");
    cplxtyPooled = pool.value<std::vector<essentia::Real>>("lowlevel.spectralComplexity");
    inharmPooled = pool.value<std::vector<essentia::Real>>("lowlevel.inharmonicity");
    barkPooled   = pool.value<std::vector<std::vector<essentia::Real>>>("lowlevel.bark");
    bfccPooled   = pool.value<std::vector<std::vector<essentia::Real>>>("lowlevel.bfcc");
    rsSpecPooled = pool.value<std::vector<std::vector<essentia::Real>>>("lowlevel.resampledSpectrum");

    
    
    std::vector<size_t> sorted_idx = sort_indexes(loudPooled);
    size_t cull_percentileLow = loudPooled.size() * percentileLow;
    std::cout << "deleting based on loudness's 20th percentile, which is " << loudPooled[sorted_idx[cull_percentileLow]] << "\n";
    sorted_idx.erase(sorted_idx.begin(), sorted_idx.begin() + cull_percentileLow);
    std::sort(sorted_idx.begin(), sorted_idx.end());
    
    for (size_t i = 0; i < sorted_idx.size(); i++){
        loudCull.push_back(loudPooled[sorted_idx[i]]);
        rmsCull.push_back(rmsPooled[sorted_idx[i]]);
        pitchCull.push_back(pitchPooled[sorted_idx[i]]);
        cplxtyCull.push_back(cplxtyPooled[sorted_idx[i]]);
        inharmCull.push_back(inharmPooled[sorted_idx[i]]);
        barkCull.push_back(barkPooled[sorted_idx[i]]);
        bfccCull.push_back(bfccPooled[sorted_idx[i]]);
        rsSpecCull.push_back(rsSpecPooled[sorted_idx[i]]);
    }
    if (loudCull.size() == rmsCull.size() == pitchCull.size() == inharmCull.size() == barkCull.size() == bfccCull.size()){
        std::cout << "sizes match\n";
    } else{
        std::cout << "SIZE MISMATCH\n";
    }
    for (size_t i = 0; i < loudCull.size(); i++){
        cullPool.add("lowlevel.loudness", loudCull[i]);
        cullPool.add("lowlevel.rms", rmsCull[i]);
        cullPool.add("lowlevel.freqs", pitchCull[i]);
        cullPool.add("lowlevel.spectralComplexity", cplxtyCull[i]);
        cullPool.add("lowlevel.inharmonicity", inharmCull[i]);
        cullPool.add("lowlevel.bark", barkCull[i]);
        cullPool.add("lowlevel.bfcc", bfccCull[i]);
        cullPool.add("lowlevel.resampledSpectrum", rsSpecCull[i]);
    }
    Algorithm* outputCull = AlgorithmFactory::create("YamlOutput",
                                                 "filename", "/Users/nicholassolem/development/audio for analysis/smooth_culled.yaml");
    outputCull->input("pool").set(cullPool);
    outputCull->compute();
    
    // aggregate the results
    essentia::Pool aggrPool, PCApool; // the pool with the aggregated values
    const char* stats[] = { "mean", "var", "min", "max" };

    Algorithm* aggr = AlgorithmFactory::create("PoolAggregator",
                                               "defaultStats", essentia::arrayToVector<std::string>(stats));

    aggr->input("input").set(pool);
    aggr->output("output").set(aggrPool);
    aggr->compute();

//    pca->input("poolIn").set(pool);
//    pca->output("poolOut").set(PCApool);
//    pca->compute();
    
    // write results to file
    std::cout << "-------- writing results to file " << outputFilename << " ---------" << std::endl;

    Algorithm* output = AlgorithmFactory::create("YamlOutput",
                                                 "filename", outputFilename);
    Algorithm* outputStats = AlgorithmFactory::create("YamlOutput",
                                                 "filename", outputFilename+"stats");
//    Algorithm* outputPCA = AlgorithmFactory::create("YamlOutput", "filename", outputFilename+"PCA");
    

    
    output->input("pool").set(pool);
    output->compute();
    
    outputStats->input("pool").set(aggrPool);
    outputStats->compute();
    
//    outputPCA->compute();

//    delete outputPCA;
//    delete pca;

    delete pyp;
    delete outputCull;
    delete outputStats;
    delete aggr;
    delete output;
//    delete rs;
//    delete zcr;
    delete rs_spec;
    delete rs;
    delete flat;
    delete inharm;
    delete hPeaks;
    delete sPeaks;
    delete sCmplxty;
    delete loud;
    delete rms;
    delete bfcc;
//    delete mfcc;
//    delete crepe;
    delete powSpec;
    delete spec;
//    delete dcBlkBig;
    delete dcBlk;
    delete w_big;
    delete w;
//    delete slicer;
    delete fc_big;
    delete fc;
    delete audio;

    essentia::shutdown();
    
    return 0;
}
