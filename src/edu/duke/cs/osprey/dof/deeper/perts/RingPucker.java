package edu.duke.cs.osprey.dof.deeper.perts;

import java.util.ArrayList;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class RingPucker extends Perturbation {

	public RingPucker(ArrayList<Residue> resDirectlyAffected) {
		super(resDirectlyAffected);
		// TODO Auto-generated constructor stub
	}

	@Override
	public boolean doPerturbationMotion(double paramVal) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Perturbation copyForNewMolecule(Molecule mol, PerturbationBlock block) {
		// TODO Auto-generated method stub
		return null;
	}

}
