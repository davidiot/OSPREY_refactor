package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.ParallelConfPartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class PartitionFunctionMinimized extends ParallelConfPartitionFunction {

	public static boolean SYNCHRONIZED_MINIMIZATION = false;
	
	protected PriorityQueue<ScoredConf> topConfs;
	protected int maxNumTopConfs;
	protected BigDecimal qstarScoreWeights;
	protected int numActiveThreads;
	protected PruningMatrix invmat;
	protected ArrayList<ScoredConf> scoredConfs;
	protected ArrayList<EnergiedConf> energiedConfs;

	public PartitionFunctionMinimized(
			EnergyMatrix emat, 
			PruningMatrix pmat, 
			PruningMatrix invmat, 
			ConfSearchFactory confSearchFactory,
			Async ecalc
			) {
		super(emat, pmat, confSearchFactory, ecalc);
		this.invmat = invmat;
		this.qstarScoreWeights = null;
		this.topConfs = null;
		this.scoredConfs = null;
		this.energiedConfs = null;
	}

	protected void writeTopConfs(int state, MSSearchProblem search) {
		if(topConfs==null || topConfs.size()==0) return;
		String seq = search.settings.getFormattedSequence();
		if(isReportingProgress) {
			System.out.println("Writing top "+ topConfs.size()+" confs:");
		}
		seq = seq.replace(" ", ".");
		String dir = "topConfs"+File.separator+"State."+state+File.separator+seq;
		ObjectIO.makeDir(dir, false);
		for(int i=topConfs.size()-1;i>-1;--i) {
			if(isReportingProgress) {
				ScoredConf head = topConfs.peek();
				double energy = head instanceof EnergiedConf ? ((EnergiedConf)head).getEnergy() : head.getScore();
				System.out.println(String.format("conf: %4d.pdf, energy: %.6f", i, energy));
			}
			String PDBFileName = dir+File.separator+i+".pdb";
			search.outputMinimizedStruct(topConfs.poll().getAssignments(), PDBFileName);
		}
	}

	protected void saveConf(ScoredConf conf) {
		if(topConfs.size() >= maxNumTopConfs) {

			ScoredConf head = topConfs.peek();
			double e1 = head instanceof EnergiedConf ? ((EnergiedConf)head).getEnergy() : head.getScore();
			double e2 = conf instanceof EnergiedConf ? ((EnergiedConf)conf).getEnergy() : conf.getScore();

			if(e1 > e2) topConfs.poll();
			else return;
		}
		topConfs.add(conf);
	}

	protected void saveEConfs(PriorityQueue<ScoredConf> other) {
		if(topConfs==null || other==null) return;
		while(other.size()>0) 
			saveConf(other.poll());
	}

	@Override
	public void init(double targetEpsilon) {

		this.targetEpsilon = targetEpsilon;

		status = Status.Estimating;
		values = new Values();

		// compute p*: boltzmann-weight the scores for all pruned conformations
		ConfSearch ptree = confSearchFactory.make(emat, invmat);
		if(ptree instanceof ConfAStarTree) ((ConfAStarTree)ptree).stopProgress();
		values.pstar = calcWeightSumUpperBound(ptree);

		// make the search tree for computing q*
		ConfSearch tree = confSearchFactory.make(emat, pmat);
		if(tree instanceof ConfAStarTree) ((ConfAStarTree)tree).stopProgress();
		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(tree);
		scoreConfs = confsSplitter.makeStream();
		energyConfs = confsSplitter.makeStream();
		numConfsEvaluated = 0;
		numConfsToScore = tree.getNumConformations();		
		qprimeUnevaluated = BigDecimal.ZERO;
		qprimeUnscored = BigDecimal.ZERO;

		qstarScoreWeights = BigDecimal.ZERO;
		numActiveThreads = 0;
		maxNumTopConfs = 0;

		stopwatch = new Stopwatch().start();
	}

	@Override
	protected BigDecimal updateQprime(EnergiedConf econf) {

		// look through the conf tree to get conf scores
		// (which should be lower bounds on the conf energy)
		while (true) {

			// read a conf from the tree
			ScoredConf conf = scoreConfs.next();
			if (conf == null) {
				qprimeUnscored = BigDecimal.ZERO;
				break;
			}

			// get the boltzmann weight
			BigDecimal scoreWeight = boltzmann.calc(conf.getScore());
			if (scoreWeight.compareTo(BigDecimal.ZERO) == 0) {
				qprimeUnscored = BigDecimal.ZERO;
				break;
			}

			// update q' parts
			numConfsToScore = numConfsToScore.subtract(BigInteger.ONE);
			qprimeUnevaluated = qprimeUnevaluated.add(scoreWeight);
			qprimeUnscored = scoreWeight.multiply(new BigDecimal(numConfsToScore));

			// stop if the bound on q' is tight enough
			double effectiveEpsilon = qprimeUnscored.divide(qprimeUnevaluated.add(qprimeUnscored), RoundingMode.HALF_UP).doubleValue();
			if (effectiveEpsilon <= 0.01) {
				break;
			}
		}

		qprimeUnevaluated = qprimeUnevaluated.subtract(boltzmann.calc(econf.getScore()));
		return qprimeUnevaluated.add(qprimeUnscored);
	}

	protected void waitForAllThreads() {
		while(numActiveThreads > 0) {
			try { this.wait(); } catch (InterruptedException e) { e.printStackTrace(); }
		}
	}

	protected ScoredConf getScoredConf(int numScored) {
		ScoredConf conf;

		if ((conf = energyConfs.next()) == null) {
			if (SYNCHRONIZED_MINIMIZATION) {
				if(status != Status.Estimated && numScored == 0) 
					status = Status.NotEnoughConformations;
			}
			else {
				waitForAllThreads();
				if(status != Status.Estimated) 
					status = Status.NotEnoughConformations;
			}
			return null;
		}

		if (boltzmann.calc(conf.getScore()).compareTo(BigDecimal.ZERO) == 0) {
			if (SYNCHRONIZED_MINIMIZATION) {
				if(status != Status.Estimated && numScored == 0) 
					status = Status.NotEnoughFiniteEnergies;
			}
			else {
				waitForAllThreads();
				if(status != Status.Estimated) 
					status = Status.NotEnoughFiniteEnergies;
			}
			return null;
		}

		++numActiveThreads;
		return conf;
	}

	protected ArrayList<ScoredConf> getScoredConfs(int numRequested) {
		ScoredConf conf;
		if(scoredConfs == null) scoredConfs = new ArrayList<>(Collections.nCopies(numRequested, null));
		if(energiedConfs == null) energiedConfs = new ArrayList<>();
		int numScored = 0;

		while(numScored < numRequested) {
			conf = getScoredConf(numScored);
			if(conf==null) break;
			scoredConfs.set(numScored++, conf);
		}

		//trim if numScored < numRequested
		if(numScored < numRequested) {
			scoredConfs.subList(numScored, scoredConfs.size()).clear();
		}

		return numScored > 0 ? scoredConfs : null;
	}

	protected void handlePhase1Conf(EnergiedConf econf) {
		if(status == Status.Estimating) {

			// get the boltzmann weight
			BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());

			// update pfunc state
			numConfsEvaluated++;
			values.qstar = values.qstar.add(energyWeight);
			values.qprime = updateQprime(econf);

			// report progress if needed
			if (isReportingProgress && numConfsEvaluated % ecalc.getParallelism() == 0) {
				phase1Output(econf);
			}

			// report confs if needed
			if (confListener != null) {
				confListener.onConf(econf);
			}

			// update status if needed
			double effectiveEpsilon = values.getEffectiveEpsilon();
			if(Double.isNaN(effectiveEpsilon)) {
				status = Status.NotEnoughFiniteEnergies;
			}
			else if (effectiveEpsilon <= targetEpsilon) {
				status = Status.Estimated;
				if (isReportingProgress) phase1Output(econf);//just to let the user know we reached epsilon
			}
		}

		--numActiveThreads;
		this.notify();
	}

	protected void handlePhase2Conf(EnergiedConf econf, BigDecimal targetScoreWeights) {
		if(status == Status.Estimating) {

			// get the boltzmann weight
			BigDecimal scoreWeight = boltzmann.calc(econf.getScore());
			qstarScoreWeights = qstarScoreWeights.add(scoreWeight);	
			BigDecimal energyWeight = boltzmann.calc(econf.getEnergy());

			// update pfunc state
			numConfsEvaluated++;
			values.qstar = values.qstar.add(energyWeight);
			values.qprime = updateQprime(econf);
			BigDecimal pdiff = targetScoreWeights.subtract(qstarScoreWeights);

			// report progress if needed
			if (isReportingProgress && numConfsEvaluated % ecalc.getParallelism() == 0) {
				phase2Output(econf, pdiff);
			}

			// report confs if needed
			if (confListener != null) {
				confListener.onConf(econf);
			}

			// update status if needed
			double effectiveEpsilon = values.getEffectiveEpsilon();
			if(Double.isNaN(effectiveEpsilon)) {
				status = Status.NotEnoughFiniteEnergies;
			}
			else if (effectiveEpsilon <= targetEpsilon) {
				status = Status.Estimated;
				if (isReportingProgress) phase2Output(econf, pdiff);
			}
		}

		--numActiveThreads;
		this.notify();
	}

	@Override
	public void compute(int maxNumConfs) {
		numActiveThreads = 0;

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		int stopAtConf = numConfsEvaluated + maxNumConfs;
		while (true) {

			// get a conf from the tree
			// lock though to keep from racing the listener thread on the conf tree
			ScoredConf conf = null;
			ArrayList<ScoredConf> scoredConfs = null;

			synchronized (this) {

				// should we keep going?
				if (!status.canContinue() || numConfsEvaluated >= stopAtConf) {
					break;
				}

				// get a list of scored confs
				if (SYNCHRONIZED_MINIMIZATION) {
					if((scoredConfs = getScoredConfs(ecalc.getParallelism())) == null)
						break;
				}

				// or a single scored conf
				else {
					if ((conf = getScoredConf(-1)) == null) 
						break;
				}
			}

			// now handle scored confs
			// wait for all confs to finish minimizing, then process confs in order of score
			if (SYNCHRONIZED_MINIMIZATION) {
				for (ScoredConf sconf : scoredConfs) ecalc.calcEnergyAsync(sconf, (EnergiedConf econf) -> {
					synchronized(this) {
						energiedConfs.add(econf);
					}
				});
				ecalc.waitForFinish();

				// sort energied confs by score
				Collections.sort(energiedConfs, new Comparator<EnergiedConf>() {
					@Override
					public int compare(EnergiedConf o1, EnergiedConf o2) {
						return o1.getScore() <= o2.getScore() ? -1 : 1;
					}
				});

				synchronized (this) {
					for (EnergiedConf econf : energiedConfs) 
						handlePhase1Conf(econf);
				}
				energiedConfs.clear();
			}

			// or do the energy calculation asynchronously
			else {
				ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {
					// energy calculation done
					// this is (potentially) running on a task executor listener thread
					// so lock to keep from racing the main thread
					synchronized (this) {
						handlePhase1Conf(econf);
					}
				});
			}
		}

		// wait for any remaining async minimizations to finish
		ecalc.waitForFinish();
	}

	void phase1Output(ScoredConf conf) {
		MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		double confVal = conf instanceof EnergiedConf ? ((EnergiedConf)conf).getEnergy() : conf.getScore();
		System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, p*: %12e, epsilon: %.6f, time: %10s, heapMem: %.0f%%",
				numConfsEvaluated, confVal, values.qstar, values.qprime, values.pstar, values.getEffectiveEpsilon(),
				stopwatch.getTime(2),
				100f*heapMem.getUsed()/heapMem.getMax()
				));
	}

	public void compute(BigDecimal targetScoreWeights, int maxNumConfs) {
		numActiveThreads = 0;

		if (!status.canContinue()) {
			throw new IllegalStateException("can't continue from status " + status);
		}

		int stopAtConf = numConfsEvaluated + maxNumConfs;
		while (true) {

			// get a conf from the tree
			// lock though to keep from racing the listener thread on the conf tree
			ScoredConf conf = null;
			ArrayList<ScoredConf> scoredConfs = null;
			synchronized (this) {

				// should we keep going?
				if (!status.canContinue() 
						|| qstarScoreWeights.compareTo(targetScoreWeights) >= 0 
						|| numConfsEvaluated >= stopAtConf) {
					break;
				}

				// get a list of scored confs
				if (SYNCHRONIZED_MINIMIZATION) {
					if((scoredConfs = getScoredConfs(ecalc.getParallelism())) == null)
						break;
				}

				// or a single scored conf
				else {
					if ((conf = getScoredConf(-1)) == null) 
						break;
				}
			}

			// wait for all confs to finish minimizing, then process confs in order of score
			if (SYNCHRONIZED_MINIMIZATION) {
				for (ScoredConf sconf : scoredConfs) ecalc.calcEnergyAsync(sconf, (EnergiedConf econf) -> {
					synchronized(this) {
						energiedConfs.add(econf);
					}
				});
				ecalc.waitForFinish();

				// sort energied confs by score
				Collections.sort(energiedConfs, new Comparator<EnergiedConf>() {
					@Override
					public int compare(EnergiedConf o1, EnergiedConf o2) {
						return o1.getScore() <= o2.getScore() ? -1 : 1;
					}
				});

				synchronized (this) {
					for (EnergiedConf econf : energiedConfs) 
						handlePhase2Conf(econf, targetScoreWeights);
				}
				energiedConfs.clear();
			}

			// or just handle the conf asynchronously
			else {
				ecalc.calcEnergyAsync(conf, (EnergiedConf econf) -> {
					// energy calculation done
					// this is (potentially) running on a task executor listener thread
					// so lock to keep from racing the main thread
					synchronized (this) {
						handlePhase2Conf(econf, targetScoreWeights);
					}
				});
			}
		}

		// wait for any remaining async minimizations to finish
		ecalc.waitForFinish();
	}

	protected void phase2Output(ScoredConf conf, BigDecimal pdiff) {
		MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		double val = conf instanceof EnergiedConf ? ((EnergiedConf)conf).getEnergy() : conf.getScore();
		System.out.println(String.format("conf: %4d, energy: %.6f, q*: %12e, q': %12e, score diff: %12e, epsilon: %.6f, time: %10s, heapMem: %.0f%%",
				numConfsEvaluated, val, values.qstar, values.qprime, pdiff, values.getEffectiveEpsilon(),
				stopwatch.getTime(2),
				100f*heapMem.getUsed()/heapMem.getMax()
				));
	}

	public void setStatus(Status val) {
		status = val;
	}

	public void cleanup() {
		scoreConfs = null;
		energyConfs = null;

		scoredConfs = null;
		energiedConfs = null;
	}
	
	public void setNumConfsEvaluated(int val) {
		numConfsEvaluated = val;
	}
	
	public int getNumConfsEvaluated() {
		return numConfsEvaluated;
	}
}
