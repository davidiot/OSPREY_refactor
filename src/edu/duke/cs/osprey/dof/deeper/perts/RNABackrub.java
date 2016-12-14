package edu.duke.cs.osprey.dof.deeper.perts;

import java.util.ArrayList;

import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.RigidBodyMotion;

public class RNABackrub extends Backrub {

	public RNABackrub(ArrayList<Residue> resDirectlyAffected) {
		super(resDirectlyAffected);
	}

	@Override
	public boolean doPerturbationMotion(double paramVal) {

		RigidBodyMotion[] rotations = calcTransRot(paramVal);
		applyModifiedBackrubLikeMotion(rotations);
		return true;// we should always be able to do an RNArub
	}

	private void applyModifiedBackrubLikeMotion(RigidBodyMotion[] motion) {
		Residue res = resDirectlyAffected.get(0);
		double[] coords = res.coords;
		// handle first plane movement (motion[0]). We move the base along with
		// the first rotation, as in backrubs
		for (String atomName : new String[] { "O4'", "C1'" }) {
			motion[0].transform(coords, res.getAtomIndexByName(atomName));
		}
		for (int atomIndex = 0; atomIndex < res.atoms.size(); atomIndex++) {
			String atomName = res.atoms.get(atomIndex).name;
			boolean isBBAtom = false;
			for (String BBAtom : HardCodedResidueInfo.possibleNABBAtoms) {
				if (atomName.equalsIgnoreCase(BBAtom)) {
					isBBAtom = true;
					break;
				}
			}
			if (!isBBAtom) {
				motion[0].transform(coords, res.getAtomIndexByName(atomName));
			}
		}
		for (String atomName : new String[] { "C2'", "H2'", "O2'", "HO2'" }) {
			motion[1].transform(coords, res.getAtomIndexByName(atomName));
		}
	}

	public double[][] extractX() {
		// We are only operating on one residue, not three
		Residue res = resDirectlyAffected.get(0);
		// x is analogous to the Calpha array in backrubs
		double[][] x = new double[3][];
		x[0] = res.getCoordsByAtomName("C3'");
		x[1] = res.getCoordsByAtomName("C1'");
		x[2] = res.getCoordsByAtomName("C4'");
		return x;
	}

	public double[] extractA() {
		Residue res = resDirectlyAffected.get(0);
		return res.getCoordsByAtomName("C2'");
	}

	public double[] extractB() {
		Residue res = resDirectlyAffected.get(0);
		return res.getCoordsByAtomName("O4'");
	}

}
