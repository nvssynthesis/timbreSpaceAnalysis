//
//  main.cpp
//  essentiaCPP_test1
//
//  Created by Nicholas Solem on 10/28/21.
//  Copyright © 2021 Nicholas Solem. All rights reserved.
//

#include <iostream>
#include "essentia.h"
#include "pool.h"
#include "algorithmfactory.h"
#include "essentiamath.h" // includes <algorithm> and <vector>
#include "noteToHertzMap.h"
//const int i=0;

using namespace essentia::standard;

int main(int argc, char* argv[]) {

//    if (argc != 3) {
//        std::cout << "ERROR: incorrect number of arguments." << std::endl;
//        std::cout << "Usage: " << argv[0] << " audio_input yaml_output" << std::endl;
//        exit(1);
//    }
    makeMap(); // build noteToHertzMap
    std::string delimiter = ".";
    std::string audioFilename = "/Users/nicholassolem/Documents/MATLAB/wtsource/uiowa/piano/Piano.pp.Db1.wav";
    std::string outputFilename = "/Users/nicholassolem/Documents/MATLAB/wtsource/uiowa/piano/Piano.pp.Db3_essentia1";
    //std::string token = audioFilename.substr(audioFilename.find(delimiter));
    //token = token.substr(token.find(delimiter));
    
    size_t pos = 0;
    std::string token[4];
    int count = 0;
    while ((pos = audioFilename.find(delimiter)) != std::string::npos) {
        token[count] = audioFilename.substr(0, pos);
        //std::cout << token[count] << std::endl;
        audioFilename.erase(0, pos + delimiter.length());
        count++;
    }
    std::string dynamic = token[1];
    std::string note = token[2];
    std::cout << noteToHertzMap[token[2]] << std::endl;

    // register the algorithms in the factory(ies)
    essentia::init();

    essentia::Pool pool;

    /////// PARAMS //////////////
    int sampleRate = 44100;
    int frameSize = 2048;
    int hopSize = 1024;
    
    AlgorithmFactory& factory = essentia::standard::AlgorithmFactory::instance();

    Algorithm* audio = factory.create("MonoLoader",
                                      "filename", audioFilename,
                                      "sampleRate", sampleRate);

    Algorithm* fc    = factory.create("FrameCutter",
                                      "frameSize", frameSize,
                                      "hopSize", hopSize);

    Algorithm* w     = factory.create("Windowing",
                                      "type", "blackmanharris62");
    Algorithm* dcBlk = factory.create("DCRemoval", "cutoffFrequency", 35); //without, there can be DC peaks which throws error for harmonic peaks
    Algorithm* spec  = factory.create("Spectrum");
    Algorithm* mfcc  = factory.create("MFCC");
    Algorithm* bfcc  = factory.create("BFCC");
    Algorithm* loud  = factory.create("Loudness");
    Algorithm* sPeaks= factory.create("SpectralPeaks");
    Algorithm* hPeaks= factory.create("HarmonicPeaks");
    Algorithm* inharm= factory.create("Inharmonicity");
    Algorithm* flat  = factory.create("Flatness");
    Algorithm* zcr   = factory.create("ZeroCrossingRate");

    /////////// CONNECTING THE ALGORITHMS ////////////////
    std::cout << "-------- connecting algos ---------" << std::endl;

    // Audio -> FrameCutter
    std::vector<essentia::Real> audioBuffer;

    audio->output("audio").set(audioBuffer);
    fc->input("signal").set(audioBuffer);   // frame creator

    // FrameCutter -> Windowing -> Spectrum
    std::vector<essentia::Real> frame, windowedFrame,hpWinFrame;
    std::vector<essentia::Real> spectrum, mfccCoeffs, mfccBands, bfccCoeffs, bfccBands,sFreqs,sMags,hFreqs,hMags;
    essentia::Real loudnesses, pitches,inharmonicities,flatnesses,zcrs;
    
    fc->output("frame").set(frame);
    w->input("frame").set(frame);

    w->output("frame").set(windowedFrame);
    //DC remove
    dcBlk->input("signal").set(windowedFrame);
    dcBlk->output("signal").set(hpWinFrame);
    spec->input("frame").set(hpWinFrame);
    spec->output("spectrum").set(spectrum);

    loud->input("signal").set(hpWinFrame);
    loud->output("loudness").set(loudnesses);

    mfcc->input("spectrum").set(spectrum);
    mfcc->output("bands").set(mfccBands);
    mfcc->output("mfcc").set(mfccCoeffs);

    bfcc->input("spectrum").set(spectrum);
    bfcc->output("bands").set(bfccBands);
    bfcc->output("bfcc").set(bfccCoeffs);

    sPeaks->input("spectrum").set(spectrum);
    sPeaks->output("frequencies").set(sFreqs);
    sPeaks->output("magnitudes").set(sMags);

    pitches = 1000.f;
    hPeaks->input("frequencies").set(sFreqs); //need to set sPeaks out first
    hPeaks->input("magnitudes").set(sMags);
    hPeaks->input("pitch").set(pitches);
    hPeaks->output("harmonicFrequencies").set(hFreqs);
    hPeaks->output("harmonicMagnitudes").set(hMags);
    
    inharm->input("frequencies").set(hFreqs);
    inharm->input("magnitudes").set(hMags);
    inharm->output("inharmonicity").set(inharmonicities);
    
    flat->input("array").set(spectrum);
    flat->output("flatness").set(flatnesses);
    
    zcr->input("signal").set(windowedFrame);
    zcr->output("zeroCrossingRate").set(zcrs);

    /////////// STARTING THE ALGORITHMS //////////////////
    std::cout << "-------- start processing " << audioFilename << " --------" << std::endl;

    audio->compute();
    
    while (true) {
        // compute a frame
        fc->compute();

        // if it was the last one (ie: it was empty), then we're done.
        if (!frame.size()) {break;}

        // if the frame is silent, just drop it and go on processing
        if (essentia::isSilent(frame)) continue;

        w->compute();
        dcBlk->compute();
        spec->compute();
        mfcc->compute();
        bfcc->compute();
        loud->compute();
        sPeaks->compute();
        hPeaks->compute();
        inharm->compute();
        flat->compute();
        zcr->compute();

        pool.add("lowlevel.mfcc", mfccCoeffs);
        pool.add("lowlevel.bfcc", bfccCoeffs);
        pool.add("lowlevel.loudness", loudnesses);
        //pool.add("lowlevel.spectralFrequencies",sFreqs);
        //pool.add("lowlevel.spectralMagnitudes",sMags);
        pool.add("lowlevel.harmonicFrequencies",hFreqs);
        pool.add("lowlevel.inharmonicity",inharmonicities);
        pool.add("lowlevel.flatness",flatnesses);
        pool.add("lowlevel.zeroCrossingRate",zcrs);
    }
    
    // aggregate the results
    essentia::Pool aggrPool; // the pool with the aggregated MFCC values
    const char* stats[] = { "mean", "var", "min", "max" };

    Algorithm* aggr = AlgorithmFactory::create("PoolAggregator",
                                               "defaultStats", essentia::arrayToVector<std::string>(stats));

    aggr->input("input").set(pool);
    aggr->output("output").set(aggrPool);
    aggr->compute();
    
    // write results to file
    std::cout << "-------- writing results to file " << outputFilename << " ---------" << std::endl;

    Algorithm* output = AlgorithmFactory::create("YamlOutput",
                                                 "filename", outputFilename);
    output->input("pool").set(pool);
    output->compute();
    
    delete audio;
    delete fc;
    delete w;
    delete spec;
    delete loud;
    delete mfcc;
    delete bfcc;
    delete sPeaks;
    delete hPeaks;
    delete inharm;
    delete flat;
    delete zcr; 
    delete output;

    essentia::shutdown();
    
    return 0;
}
