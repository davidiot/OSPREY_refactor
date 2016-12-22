package edu.duke.cs.osprey.dof.deeper.perts;

import java.util.ArrayList;
import java.util.HashSet;

import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.VectorAlgebra;

/**
 * A perturbation for RNA residues that allows the base and OH groups to move
 * 
 * @author dzhou
 */
public class RingPucker extends Perturbation {

	// For O(1) lookup when checking which atoms to move
	private static HashSet<String> atomsToHoldInPlace;

	// initialize the HashSet
	static {
		atomsToHoldInPlace = new HashSet<String>();
		for (String atom : HardCodedResidueInfo.possibleNABBAtoms) {
			if (!atom.equals("C1'")) {
				atomsToHoldInPlace.add(atom);
			}
		}
	}

	// Constructor
	public RingPucker(ArrayList<Residue> resDirectlyAffected) {
		super(resDirectlyAffected);
	}

	@Override
	public boolean doPerturbationMotion(double paramVal) {
		// Apply the perturbation
		// Use an arbitrary param (primary ring pucker angle in degrees)
		// Don't store rotation matrices or translations
		Residue res = resDirectlyAffected.get(0);
		double originalDihedralAngle = Protractor.measureDihedral(res.coords, res.getAtomIndexByName("C1'"),
				res.getAtomIndexByName("C2'"), res.getAtomIndexByName("C3'"), res.getAtomIndexByName("C4'"));
		rotateC1(res, paramVal);
		adjustC2(res, originalDihedralAngle);

		return true; // we should always be able to do a ring pucker in RNA
	}

	/**
	 * Rotate C1 (the carbon connected to the nitrogenous base), as well as the
	 * base and H1'. We rotate about the axis defined by C2' and O4' by the
	 * angle specified by paramVal
	 */
	private void rotateC1(Residue res, double paramVal) {
		double[] rotationAxis = VectorAlgebra.subtract(res.getCoordsByAtomName("C2'"), res.getCoordsByAtomName("O4'"));
		RigidBodyMotion rbm = new RigidBodyMotion(res.getCoordsByAtomName("C2'"), rotationAxis, paramVal, false);

		for (int atomIndex = 0; atomIndex < res.atoms.size(); atomIndex++) {
			String atomName = res.atoms.get(atomIndex).name;
			if (!atomsToHoldInPlace.contains(atomName)) {
				rbm.transform(res.coords, atomIndex);
			}
		}

	}

	/**
	 * Adjust C2 (the carbon connected to the OH group) so that it is close to
	 * tetahedral again.
	 */
	private void adjustC2(Residue res, double originalDihedralAngle) {
		double newDihedralAngle = Protractor.measureDihedral(res.coords, res.getAtomIndexByName("C1'"),
				res.getAtomIndexByName("C2'"), res.getAtomIndexByName("C3'"), res.getAtomIndexByName("C4'"));

		double compensationAngle = originalDihedralAngle - newDihedralAngle;

		double[] rotationAxis = VectorAlgebra.subtract(res.getCoordsByAtomName("C3'"), res.getCoordsByAtomName("C2'"));
		RigidBodyMotion rbm = new RigidBodyMotion(res.getCoordsByAtomName("C3'"), rotationAxis, compensationAngle,
				false);

		for (String atomName : new String[] { "H2'", "O2'", "HO2'" }) {
			rbm.transform(res.coords, res.getAtomIndexByName(atomName));
		}

	}

	@Override
	public Perturbation copyForNewMolecule(Molecule mol, PerturbationBlock block) {
		RingPucker rp = new RingPucker(Residue.equivalentInMolec(resDirectlyAffected, mol));
		rp.curParamVal = curParamVal;
		rp.indexInBlock = indexInBlock;
		rp.block = block;
		return rp;
	}

}
