package edu.duke.cs.osprey.control;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.multistatekstar.InputValidation;
import edu.duke.cs.osprey.multistatekstar.KStarScore;
import edu.duke.cs.osprey.multistatekstar.LMB;
import edu.duke.cs.osprey.multistatekstar.MSConfigFileParser;
import edu.duke.cs.osprey.multistatekstar.MSKStarFactory;
import edu.duke.cs.osprey.multistatekstar.MSKStarTree;
import edu.duke.cs.osprey.multistatekstar.MSSearchProblem;
import edu.duke.cs.osprey.multistatekstar.MSSearchSettings;
import edu.duke.cs.osprey.multistatekstar.PartitionFunctionMinimized;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.StringParsing;

public class MSKStarDoer {

	MSKStarTree tree;//tree used for the MultiStateKStar search
	int numSeqsWanted;//number of desired sequences
	int numMaxMut;//max number of mutations from wt

	LMB objFcn;//objective function for MultiStateKStar search, i.e. the f-score
	LMB[] msConstr;//global constraints for search. (partial)sequences that violate constraints are pruned
	LMB[][] sConstr;//state-specific constraints
	int numStates;//number of states considered
	int numMutRes;//number of mutable positions

	ArrayList<String[]> wtSeqs;//bound state wild type sequences for each state

	ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions;//AA types allowed at 
	//each mutable position for each substate

	ArrayList<ArrayList<ArrayList<Integer>>> state2MutableResNums;
	//For each state, a list of flexible/mutable positions
	//these will be listed directly in system cfg files

	String stateArgs[][];//search arguments for each state

	ParamSet msParams;//multistate spec params
	MSConfigFileParser[] cfps;//config file parsers for each state

	SearchProblem[][] searchDisc;//continuous search problems
	SearchProblem[][] searchCont;//discrete search problems

	ConfEnergyCalculator.Async[][] ecalcsCont;//global continuous energy calculator objects
	ConfEnergyCalculator.Async[][] ecalcsDisc;//global discrete energy calculator objects

	public MSKStarDoer(String args[]) {

		// silence warnings when using non-amino acids
		BigForcefieldEnergy.ParamInfo.printWarnings = false;
		ForcefieldParams.printWarnings = false;

		//fill in all the settings
		//each state will have its own config file parser

		System.out.println();
		System.out.println("Performing multistate K*");
		System.out.println();

		//check format of args
		if(!args[0].equalsIgnoreCase("-c"))
			throw new RuntimeException("ERROR: bad arguments (should start with -c)");

		// multistate spec parameters
		msParams = new ParamSet();
		msParams.setVerbosity(false);
		msParams.addParamsFromFile(args[4]);//read multistate parameters
		msParams.addDefaultParams();

		numSeqsWanted = msParams.getInt("NUMSEQS");
		numStates = msParams.getInt("NUMSTATES");
		numMutRes = msParams.getInt("NUMMUTRES");
		numMaxMut = msParams.getInt("NUMMAXMUT");
		int numConstr = msParams.getInt("NUMSTATECONSTR");

		stateArgs = new String[numStates][];

		objFcn = new LMB(msParams.getValue("OBJFCN"), numStates);

		// might need to adjust these with wt later
		msConstr = new LMB[numConstr];
		for(int constr=0; constr<numConstr; constr++)
			msConstr[constr] = new LMB(msParams.getValue("STATECONSTR"+constr), numStates);

		sConstr = new LMB[numStates][];

		cfps = new MSConfigFileParser[numStates];

		searchCont = new SearchProblem[numStates][];
		searchDisc = new SearchProblem[numStates][];

		ecalcsCont = new ConfEnergyCalculator.Async[numStates][];
		ecalcsDisc = new ConfEnergyCalculator.Async[numStates][];

		System.out.println();
		System.out.println("Checking multistate K* parameters for consistency");
		System.out.println();

		state2MutableResNums = new ArrayList<>();
		AATypeOptions = new ArrayList<>();
		wtSeqs = new ArrayList<>();

		InputValidation inputValidation = new InputValidation(AATypeOptions, state2MutableResNums);
		inputValidation.handleObjFcn(msParams, objFcn);
		inputValidation.handleConstraints(msParams, msConstr);

		for(int state=0; state<numStates; state++) {

			System.out.println();
			System.out.println("Checking state "+state+" parameters");
			System.out.println();

			cfps[state] = makeStateCfp(state);
			cfps[state].getParams().setValue("PRUNEPARTIALSEQCONFS", msParams.getValue("PRUNEPARTIALSEQCONFS"));
			inputValidation.handleDoMinimize(cfps);

			ParamSet sParams = cfps[state].getParams();
			inputValidation.handleStateParams(state, sParams, msParams);
			state2MutableResNums.add(stateMutableRes(state, cfps[state], numMutRes));

			for(int subState=0; subState<state2MutableResNums.get(state).size(); ++subState){
				inputValidation.handleAATypeOptions(state, subState, cfps[state]);
				//get bound substate wt sequence
				if(subState==state2MutableResNums.get(state).size()-1)
					wtSeqs.add(cfps[state].getWtSeq(state2MutableResNums.get(state).get(subState)));
			}

			//populate state-specific constraints
			int numUbConstr = sParams.getInt("NUMSTRANDCONSTR");
			int numPartFuncs = sParams.getInt("NUMOFSTRANDS")+1;
			sConstr[state] = new LMB[numUbConstr];
			for(int constr=0;constr<numUbConstr;constr++)
				sConstr[state][constr] = new LMB(sParams.getValue("STRANDCONSTR"+constr), numPartFuncs);

			System.out.println();
			System.out.println("State "+state+" parameters checked");
			System.out.println();
		}

		state2MutableResNums.trimToSize();
		wtSeqs.trimToSize();

		System.out.println();
		System.out.println("Preparing search problems and matrices for multistate K*");
		System.out.println();

		for(int state=0; state<numStates; state++) {

			searchCont[state] = makeStateSearchProblems(state, true, cfps[state]);//continuous flex
			searchDisc[state] = makeStateSearchProblems(state, false, cfps[state]);//discrete flex

			//set verbosity to false now
			cfps[state].getParams().setVerbosity(false);
			
			System.out.println();
			System.out.println("State "+state+" matrices ready");
			System.out.println();
		}
	}

	protected ConfEnergyCalculator.Async[][] makeEnergyCalculators(boolean cont) {
		ConfEnergyCalculator.Async[][] ans = new ConfEnergyCalculator.Async[numStates][];
		for(int state=0;state<numStates;++state) {
			ans[state] = makeEnergyCalculators(state, cont);
		}
		return ans;
	}

	private ConfEnergyCalculator.Async[] makeEnergyCalculators(int state, boolean cont) {
		SearchProblem[] search = cont ? searchCont[state] : searchDisc[state];
		Parallelism parallelism = cont ? Parallelism.makeFromConfig(cfps[state]) : Parallelism.makeDefault();
		ConfEnergyCalculator.Async[] ans = new ConfEnergyCalculator.Async[search.length];
		for(int substate=0;substate<search.length;++substate) {
			ans[substate] = MSKStarFactory.makeEnergyCalculator(cfps[state], search[substate], parallelism);
		}
		return ans;
	}

	private void cleanupEnergyCalculators(ConfEnergyCalculator.Async[][] ecalcs, int state) {
		if(ecalcs[state]==null) return;
		for(int substate=0;substate<ecalcs[state].length;++substate) {
			ConfEnergyCalculator.Async ecalc = ecalcs[state][substate];
			if(ecalc != null) {
				ecalcs[state][substate].cleanup();
				ecalcs[state][substate] = null;
			}
		}
		ecalcs[state] = null;
	}

	private void cleanupEnergyCalculators(ConfEnergyCalculator.Async[][] ecalcs) {
		if(ecalcs == null) return;
		for(int state=0;state<ecalcs.length;++state) cleanupEnergyCalculators(ecalcs, state);
		ecalcs = null;
	}

	private void cleanup() {
		cleanupEnergyCalculators(ecalcsCont);
		cleanupEnergyCalculators(ecalcsDisc);
	}

	private String getStateKStarScore(int state, ArrayList<String> boundStateAATypes) {

		//get arraylist formatted sequence for each substate
		ArrayList<ArrayList<ArrayList<String>>> subStateAATypes = new ArrayList<>();
		for(int subState=0;subState<state2MutableResNums.get(state).size();++subState){
			//state2MutableResNums.get(state).get(subState)) contains substate flexible residues
			//last substate is the bound state
			subStateAATypes.add(new ArrayList<>());
			ArrayList<Integer> subStateResNums = state2MutableResNums.get(state).get(subState);
			ArrayList<Integer> boundStateResNums = state2MutableResNums.get(state).get(state2MutableResNums.get(state).size()-1);
			for(int resNum : subStateResNums){
				int index = boundStateResNums.indexOf(resNum);
				ArrayList<String> aa = new ArrayList<>(); aa.add(boundStateAATypes.get(index));
				subStateAATypes.get(subState).add(aa);
			}
		}

		ParamSet sParams = cfps[state].getParams();
		int numPartFuncs = sParams.getInt("NUMOFSTRANDS")+1;
		boolean doMinimize = sParams.getBool("DOMINIMIZE");
		KStarScoreType scoreType = MSKStarFactory.getKStarScoreType(sParams);

		//populate search problems
		MSSearchProblem[] singleSeqSearchCont = doMinimize ? new MSSearchProblem[numPartFuncs] : null;
		MSSearchProblem[] singleSeqSearchDisc = !doMinimize ? new MSSearchProblem[numPartFuncs] : null;
		for(int subState=0;subState<numPartFuncs;++subState){

			MSSearchSettings spSet = new MSSearchSettings();
			spSet.AATypeOptions = subStateAATypes.get(subState);
			ArrayList<String> mutRes = new ArrayList<>();
			for(int i:state2MutableResNums.get(state).get(subState)) mutRes.add(String.valueOf(i));
			spSet.mutRes = mutRes;
			spSet.stericThreshold = sParams.getDouble("STERICTHRESH");
			spSet.pruningWindow = sParams.getDouble("IVAL") + sParams.getDouble("EW");

			if(doMinimize)
				singleSeqSearchCont[subState] = new MSSearchProblem(searchCont[state][subState], spSet);
			else
				singleSeqSearchDisc[subState] = new MSSearchProblem(searchDisc[state][subState], spSet);
		}

		KStarScore score = MSKStarFactory.makeKStarScore(
				msParams, state, cfps[state], sConstr[state],
				singleSeqSearchCont, singleSeqSearchDisc,
				ecalcsCont[state], ecalcsDisc[state], scoreType
				);
		score.compute(Integer.MAX_VALUE);
		return score.toString();
	}

	/**
	 * creates and prunes all-seq energy matrices for each substate within a state
	 * @param state
	 * @param cont
	 * @param stateCfp
	 * @return
	 */
	private SearchProblem[] makeStateSearchProblems(int state, boolean cont, 
			MSConfigFileParser stateCfp) {

		ParamSet sParams = stateCfp.getParams();
		int numUbStates = sParams.getInt("NUMOFSTRANDS");
		String flexibility = cont ? "continuous" : "discrete";

		SearchProblem[] subStateSps = new SearchProblem[numUbStates+1];

		for(int subState=0;subState<subStateSps.length;++subState) {
			subStateSps[subState] = stateCfp.getSearchProblem(state, subState, 
					state2MutableResNums.get(state).get(subState), cont);

			//make emats
			subStateSps[subState].loadEnergyMatrix();

			//prune
			if(!sParams.getBool("UsePoissonBoltzmann")) {
				PruningControl pc = stateCfp.setupPruning(subStateSps[subState], 
						sParams.getDouble("Ival")+sParams.getDouble("Ew"), 
						sParams.getBool("UseEpic"), 
						sParams.getBool("UseTupExp"));
				//silence output
				pc.setReportMode(null);
				pc.prune();
			}

			System.out.println();
			System.out.println("State "+state+"."+subState+" "+flexibility+" matrix ready");
			System.out.println();
		}

		return subStateSps;
	}

	/**
	 * creates mutable residue list for all substates in all states
	 * @param state
	 * @param stateCfp
	 * @param numTreeLevels
	 * @return
	 */
	private ArrayList<ArrayList<Integer>> stateMutableRes(int state, MSConfigFileParser stateCfp, int numTreeLevels){
		ParamSet sParams = stateCfp.getParams();
		int numUbStates = sParams.getInt("NUMOFSTRANDS");
		ArrayList<ArrayList<Integer>> m2s = new ArrayList<>();
		for(int ubState=0;ubState<=numUbStates;++ubState) m2s.add(new ArrayList<>());

		for(int ubState=0;ubState<numUbStates;++ubState) {
			StringTokenizer st = new StringTokenizer(sParams.getValue("STRANDMUT"+ubState));
			while(st.hasMoreTokens()) m2s.get(ubState).add(Integer.valueOf(st.nextToken()));
			//append to complex residues
			m2s.get(numUbStates).addAll(m2s.get(ubState));
			m2s.get(numUbStates).trimToSize();
		}

		if(m2s.get(numUbStates).size()!=numTreeLevels){
			throw new RuntimeException("ERROR: SeqTree has "+numTreeLevels+" mutable positions "
					+" but "+m2s.size()+" are listed for state "+state);
		}

		m2s.trimToSize();
		return m2s;
	}

	/**
	 * makes a config file parser for each state
	 * @param state
	 * @return
	 */
	private MSConfigFileParser makeStateCfp(int state) {
		//We expect input of the form KStar0.cfg System0.cfg DEE0.cfg
		String stateConfigFiles = msParams.getValue("STATECFGFILES"+state);
		String stateKStFile = StringParsing.getToken(stateConfigFiles, 1);
		String stateSysFile = StringParsing.getToken(stateConfigFiles, 2);
		String stateDEEFile = StringParsing.getToken(stateConfigFiles, 3);

		stateArgs[state] = new String[] {"-c", stateKStFile, "n/a", stateSysFile, stateDEEFile};

		MSConfigFileParser stateCfp = new MSConfigFileParser(stateArgs[state]);
		//verbosity is false by default; too many arguments clutter the display
		stateCfp.loadData();

		return stateCfp;
	}

	protected void printAllSeqs(){
		ArrayList<ArrayList<ArrayList<String>>> stateSeqLists = listAllSeqs();
		for(int state=0;state<stateSeqLists.size();++state){

			int numSeqs=stateSeqLists.get(state).size();

			System.out.println();
			System.out.println("State"+state+": "+numSeqs+" sequences with <= "+numMaxMut+" mutation(s) from wild-type");
			System.out.println();

			for(ArrayList<String> seq : stateSeqLists.get(state)){ 
				for(String aa : seq) System.out.print(aa+" ");
				System.out.println();
			}
		}
	}

	/**
	 * Get a list of all possible sequences that are maxNumMut mutations away
	 * from the wilt type sequence.
	 * @return
	 */
	private ArrayList<ArrayList<ArrayList<String>>> listAllSeqs(){

		System.out.println();
		System.out.print("Counting number of possible sequences...");

		//for all bound states, list all possible sequences for the mutable residues,
		//based on AATypeOptions
		ArrayList<ArrayList<ArrayList<String>>> ans = new ArrayList<>();

		//pre-allocate buffer
		String[] buf = new String[numMutRes];

		for(int state=0;state<numStates;++state) {
			int subState = AATypeOptions.get(state).size()-1;
			ArrayList<ArrayList<String>> subStateAATypeOptions = AATypeOptions.get(state).get(subState);
			ArrayList<ArrayList<String>> stateOutput = new ArrayList<>();

			//get allowed sequences for this state's bound complex
			listAllSeqsHelper(subStateAATypeOptions, stateOutput, wtSeqs.get(state), buf, 0, 0);
			stateOutput.trimToSize();
			ans.add(stateOutput);
		}

		System.out.println("done");
		System.out.println("Number of possible sequences: "+ans.get(0).size());
		System.out.println();

		ans.trimToSize();
		return ans;
	}

	private void listAllSeqsHelper(ArrayList<ArrayList<String>> subStateAATypeOptions, 
			ArrayList<ArrayList<String>> stateOutput, String[] wt, String[] buf, int depth, int dist){
		//List all sequences for the subset of mutable positions with max distance
		//from wt starting at depth=0 and going to the last mutable position
		if(depth==numMutRes){
			//String[] seq = new String[numTreeLevels];
			//System.arraycopy(buf, 0, seq, 0, numTreeLevels);
			ArrayList<String> seq = new ArrayList<String>(Arrays.asList(buf));
			seq.trimToSize();
			stateOutput.add(seq);
			return;
		}

		for(int aaIndex=0;aaIndex<subStateAATypeOptions.get(depth).size();++aaIndex){
			buf[depth]=subStateAATypeOptions.get(depth).get(aaIndex);
			int nDist=buf[depth].equalsIgnoreCase(wt[depth]) ? dist : dist+1;
			if(nDist>numMaxMut) continue;
			listAllSeqsHelper(subStateAATypeOptions, stateOutput, wt, buf, depth+1, nDist);
		}
	}

	public void calcBestSequences() {
		PartitionFunctionMinimized.SYNCHRONIZED_MINIMIZATION = this.msParams.getBool("SYNCHRONIZEDMINIMIZATION");
		final String algOption = msParams.getValue("MultStateAlgOption");
		switch(algOption.toLowerCase()) {
		case "exhaustive":
			exhaustiveMultistateSearch(msParams.getDouble("TIMEOUT"));
			return;
		case "sublinear":
			subLinearMultiStateSearch();
			return;
		default:
			throw new UnsupportedOperationException("ERROR: "+algOption+" is not supported for MULTISTATEALGOPTION");
		}
	}

	public void subLinearMultiStateSearch() {

		System.out.println();
		System.out.println("Performing sub-linear multistate K*");
		System.out.println();

		//create energy calculators
		boolean doMinimize = false;
		for(int state=0;state<numStates;++state) {
			doMinimize = cfps[state].getParams().getBool("DOMINIMIZE");
			if(doMinimize) ecalcsCont[state] = makeEnergyCalculators(state, true);
			ecalcsDisc[state] = makeEnergyCalculators(state, false);
		}

		tree = new MSKStarTree(numMutRes, numStates, numMaxMut, numSeqsWanted, objFcn, 
				msConstr, sConstr, state2MutableResNums, AATypeOptions, wtSeqs, 
				searchCont, searchDisc, ecalcsCont, ecalcsDisc, msParams, cfps);

		ArrayList<String> bestSequences = new ArrayList<>();
		String fname = "sequences-sublinear."+msParams.getValue("RUNNAME")+".txt";
		ObjectIO.delete(fname);
		PrintStream fout = null;
		String elapsed = null;
		Stopwatch stopwatch = new Stopwatch().start();

		try {
			fout = new PrintStream(new FileOutputStream(new File(fname), true));

			for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
				String seq = tree.nextSeq();//this will find the best sequence
				if(seq == null)//empty sequence...indicates no more sequence possibilities
					break;
				else {
					bestSequences.add(seq);
					fout.println(seq); fout.flush();
				}
			}

			elapsed = stopwatch.getTime(2);

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if(fout != null) {fout.flush();fout.close();}
			cleanup();//free energy calculators
		}

		System.out.println();
		System.out.println("Finished sub-linear multi-state K* in "+elapsed);
		System.out.println();

		//this prints out the total number of sequences
		listAllSeqs();
	}

	/**
	 * Verify algorithm results by doing exhaustive search
	 */
	private void exhaustiveMultistateSearch(double timeoutHrs) {

		System.out.println();
		System.out.println("Checking multi-state K* by exhaustive search");
		System.out.println();

		Stopwatch stopwatch = new Stopwatch().start();
		ArrayList<ArrayList<ArrayList<String>>> seqList = listAllSeqs();

		//process selected mutants
		String mutFname = msParams.getValue("MUTFILE", "");
		if(mutFname.length()>0) {
			//extract mutations from file
			ArrayList<ArrayList<String>> mutSeqs = readMutFile(mutFname);
			//make sure they exist in seqList; print warnings for/remove those that don't
			ArrayList<ArrayList<String>> keep = new ArrayList<>();
			for(ArrayList<String> seq : mutSeqs){
				if(seqList.get(0).contains(seq)) keep.add(seq);
			}
			keep.trimToSize();
			mutSeqs.removeAll(keep);
			if(mutSeqs.size()>0){
				System.out.println();
				System.out.println("WARNING: the following "+ mutSeqs.size() +" sequence(s) in MUTFILE "
						+ "were not in the original sequence space and will be ignored");
				for(ArrayList<String> seq : mutSeqs){
					for(String str : seq) System.out.print(str + " ");
					System.out.println();
				}
				System.out.println();
			}
			//set seqlist to the mutants of interest
			for(int state=0;state<numStates;++state) seqList.set(state, keep);
		}

		//resume
		String fname = "sequences-exhaustive."+msParams.getValue("RUNNAME")+".txt";
		boolean resume = msParams.getBool("RESUME");
		if(resume) {
			ArrayList<ArrayList<ArrayList<String>>> completed = getCompletedSeqs(fname);
			for(int state=0;state<numStates;++state){
				seqList.get(state).removeAll(completed.get(state));
			}
		}

		String[][] stateKSS = new String[numStates][];
		for(int state=0;state<numStates;++state) stateKSS[state] = new String[seqList.get(state).size()];
		PrintStream fout = null;
		try {
			if(!resume) 
				ObjectIO.delete(fname);

			fout = new PrintStream(new FileOutputStream(new File(fname), true));

			for(int state=0; state<numStates; state++){
				boolean doMinimize = cfps[state].getParams().getBool("DOMINIMIZE");

				if(stateKSS[state].length>0) {
					fout.println();
					fout.println("State"+state+": ");
					fout.println();

					System.out.println();
					System.out.println("State"+state+": ");
					System.out.println();
				}

				//make energy calculators for this state
				if(doMinimize) ecalcsCont[state] = makeEnergyCalculators(state, true);
				else ecalcsDisc[state] = makeEnergyCalculators(state, false);

				for(int seqNum=0; seqNum<stateKSS[state].length; seqNum++){
					System.out.println();
					System.out.print("Computing sequence "+(seqNum+1)+"/"+stateKSS[state].length+"...");
					stateKSS[state][seqNum] = getStateKStarScore(state, seqList.get(state).get(seqNum));
					fout.println(stateKSS[state][seqNum]); fout.flush();
					System.out.println("done, elapsed: "+stopwatch.getTime(2));
					System.out.println();
					
					if(timeoutHrs > 0 && stopwatch.getTimeH() >= timeoutHrs) {
						System.out.println("WARNING: running time has exceeded "
								+ "the allotted "+timeoutHrs+" hours. Stopping execution.");
						break;
					}
				}

				searchCont[state] = null;
				searchDisc[state] = null;
				cleanupEnergyCalculators(ecalcsCont, state);
				cleanupEnergyCalculators(ecalcsDisc, state);
			}

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if(fout!=null){fout.flush(); fout.close();}
			cleanup();
		}

		printAllKStarScores(stateKSS);

		System.out.println();
		System.out.println("Finished checking multi-state K* by exhaustive search in "+stopwatch.getTime(2));
		System.out.println();
	}

	/**
	 * Prints all K* scores for all sequences in all states
	 * @param stateKSS
	 */
	private void printAllKStarScores(String[][] stateKSS){
		for(int state=0;state<numStates;++state){
			if(stateKSS[state].length>0) {
				System.out.println();
				System.out.println("State"+state+": ");
				System.out.println();
			}
			String[] kss = stateKSS[state];
			for(int subState=0;subState<kss.length;++subState){
				System.out.println(kss[subState]);
			}
		}
	}

	private ArrayList<ArrayList<String>> readMutFile(String fname) {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();

		if(!(new File(fname)).exists()) return ans;

		try (BufferedReader br = new BufferedReader(new FileReader(fname))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if(line.length()==0) continue;

				ArrayList<String> seq = new ArrayList<>();
				StringTokenizer st = new StringTokenizer(line);
				while(st.hasMoreTokens()) seq.add(st.nextToken());
				seq.trimToSize();
				ans.add(seq);
			}
		} catch (IOException e) { e.printStackTrace(); }

		ans.trimToSize();
		return ans;
	}

	private ArrayList<ArrayList<ArrayList<String>>> getCompletedSeqs(String fname) {
		ArrayList<ArrayList<ArrayList<String>>> ans = new ArrayList<>();
		for(int state=0;state<numStates;++state) ans.add(new ArrayList<>());

		if(!(new File(fname)).exists()) return ans;

		for(int state=0;state<numStates;++state) ans.add(new ArrayList<>());

		try (BufferedReader br = new BufferedReader(new FileReader(fname))) {
			String line;
			int state=0;
			while ((line = br.readLine()) != null) {

				if(line.length()==0) continue;

				//first find state
				else if(line.startsWith("State")) {
					line = line.toLowerCase();
					line = line.replace("state", "");
					line = line.replace(":", "").trim();
					state = Integer.valueOf(line);
					continue;
				}

				//records for state
				else if(line.startsWith("Seq")) {
					StringTokenizer st = new StringTokenizer(line, ":");
					while(st.hasMoreTokens()) {
						String token = st.nextToken();
						if(!token.contains("score")) continue;
						token = token.replace("score", "");
						token = token.replace(",", "");
						token = token.trim();

						StringTokenizer st1 = new StringTokenizer(token);
						ArrayList<String> val1 = new ArrayList<>();
						while(st1.hasMoreTokens()) val1.add(st1.nextToken().split("-")[0]);
						val1.trimToSize();
						ans.get(state).add(val1);
					}
				}
			}

			br.close();

		} catch (IOException e) { e.printStackTrace(); }

		return ans;
	}
}
