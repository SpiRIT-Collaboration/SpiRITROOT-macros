void GetParameters(Int_t run, Int_t &numTotal, TString &GCData, TString &GGData)
{
  Int_t r, n;
  ifstream file("run_info.txt");
  while (file >> r >> n) {
    if (r == run) {
      numTotal = n; 
      GCData = "";
      Int_t ggRun = 0;

      // Sn108
      if (run >= 2235 && run <= 2238) ggRun = 2229; // Don't use them, both runs and GG noise data. Peaking time 117ns.
      if (run >= 2240 && run <= 2258) ggRun = 2244; // Peaking time 232ns. Probably too much work to correct those.
      if (run >= 2259 && run <= 2289) ggRun = 2263; // Peaking time 117ns for the runs below.
      if (run >= 2290 && run <= 2320) ggRun = 2290;
      if (run >= 2321 && run <= 2434) ggRun = 2321; // Justin said 2321 is the best since there's lots of problem in GG at that time.
      if (run >= 2435 && run <= 2453) ggRun = 2435;
      if (run >= 2461 && run <= 2484) ggRun = 2464;
      if (run >= 2498 && run <= 2509) ggRun = 2504;

      // Sn112+Sn124 - Not confirmed. Just based on the reconstruction sheet.
      if (run >= 2520 && run <= 2575) ggRun = 2545;
      if (run >= 2577 && run <= 2647) ggRun = 2576;
      if (run >= 2649 && run <= 2653) ggRun = 2648;

      // Sn124+Sn112
      if (run >= 3058 && run <= 3085) ggRun = 3060;
      if (run >= 3086 && run <= 3095) ggRun = 3086;
      if (run >= 3096 && run <= 3100) ggRun = 3096;
      if (run >= 3101 && run <= 3103) ggRun = 3101;
      if (run >= 3104 && run <= 3130) ggRun = 3104;
      if (run >= 3131 && run <= 3146) ggRun = 3131;
      if (run >= 3147 && run <= 3162) ggRun = 3147;
      if (run >= 3163 && run <= 3177) ggRun = 3163;
      if (run >= 3178 && run <= 3184) ggRun = 3178;

      // Sn132
      if (run >= 2841 && run <= 2852) ggRun = 2842;
      if (run >= 2855 && run <= 2861) ggRun = 2853;
      if (run >= 2877 && run <= 2884) ggRun = 2876;
      if (run >= 2887 && run <= 2894) ggRun = 2886;
      if (run >= 2896 && run <= 2905) ggRun = 2895;
      if (run >= 2907 && run <= 2917) ggRun = 2906;
      if (run >= 2919 && run <= 2927) ggRun = 2918;
      if (run >= 2929 && run <= 2937) ggRun = 2928;
      if (run >= 2939 && run <= 2946) ggRun = 2938;
      if (run >= 2948 && run <= 2962) ggRun = 2947;
      if (run >= 2964 && run <= 2973) ggRun = 2963;
      if (run >= 2975 && run <= 2986) ggRun = 2974;
      if (run >= 2988 && run <= 2997) ggRun = 2987;
      if (run >= 2999 && run <= 3010) ggRun = 2998;
      if (run >= 3039 && run <= 3039) ggRun = 3037;

      // Modify the path of pulser calibration file. It's in parameters folder in the latest SpiRITROOT.
      GCData = "/mnt/spirit/rawdata/misc/gainCalibration_groundPlane_120fC_117ns_20160509.root"; // Fix me

      // Change the gating grid file path
      if (ggRun != 0)
        GGData = Form("/mnt/spirit/rawdata/misc/ggNoise/ggNoise_%d.root",ggRun); // Fix me
      else 
        GGData = "";
      break; 
    }
  }
  file.close();
}

void create_submit(Int_t run = 2900)
{ 
  Int_t numSplit = 2000; // Fix me if you want to change the event number in a split
  Int_t numTotal = 0;
  TString GCData, GGData;
  GetParameters(run, numTotal, GCData, GGData);
  if (numTotal == 0) 
    continue;

  Int_t mjs = (numTotal/3+1)/numSplit;

  TString version = ".";
  TString fileName = Form("submit/r%d",run) + version + "sh";
  TString spiritrootDIR = "/mnt/spirit/analysis/changj/SpiRITROOT.latest"; // Fix me
  ofstream out(fileName);

  // This preamble is for ember.compute. This should be modified for your system. Fix me
  out << "#!/usr/bin/env bash" << endl;
  out << "#--- sbatch option ---#" << endl;
  out << "#SBATCH --ntasks=1" << endl;
  out << "#SBATCH --cpus-per-task=8" << endl;
  out << "#SBATCH --mem-per-cpu=2000" << endl;
  out << Form("#SBATCH --array=0-%d", mjs) << endl; 
  out << endl;
  // Up to here

  out << "source " + spiritrootDir + "/build/config.sh" << endl;
  out << "cd " + spiritrootDir + "/macros/" << endl; 
  out << endl;
  out << Form("RUN=%d", run) << endl;
  out << Form("NTOTAL=%d", numTotal) << endl;
  out << Form("NSPLIT=%d", numSplit) << endl; 
  out << "GCData=" << GCData << endl;
  out << "GGData=" << GGData << endl;
  out << endl;
  // You must replace "SLURM_ARRAY_TASK" into some variable in your system. Fix me
  out << "SPLIT=$((3*SLURM_ARRAY_TASK_ID+0)); root run_reco_experiment.C\\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\\\"$GCData\\\",\\\"$GGData\\\"\\) -b -q -l > log/log_run$RUN\\_$SPLIT.log 2>&1 &" << endl;
  out << "SPLIT=$((3*SLURM_ARRAY_TASK_ID+1)); root run_reco_experiment.C\\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\\\"$GCData\\\",\\\"$GGData\\\"\\) -b -q -l > log/log_run$RUN\\_$SPLIT.log 2>&1 &" << endl;
  out << "SPLIT=$((3*SLURM_ARRAY_TASK_ID+2)); root run_reco_experiment.C\\($RUN,$NTOTAL,$SPLIT,$NSPLIT,\\\"$GCData\\\",\\\"$GGData\\\"\\) -b -q -l > log/log_run$RUN\\_$SPLIT.log 2>&1 &" << endl;
  out << endl;
  out << "wait" << endl;

  cout << fileName << " " << numTotal << " " << GCData << " " << GGData << endl;
  out.close();
}
