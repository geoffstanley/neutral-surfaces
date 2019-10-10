/*
 *	Copyright (C) 2012 Visualization & Graphics Lab (VGL), Indian Institute of Science
 *
 *	This file is part of Recon, a library to compute Reeb graphs.
 *
 *	Recon is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	Recon is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public License
 *	along with Recon.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Author(s):	Harish Doraiswamy
 *	Version	 :	1.0
 *
 *	Modified by : -- 
 *	Date : --
 *	Changes  : --
 */
package vgl.iisc.recon.incore;



public class DisjointSetsInt {

	int [] set;

	public DisjointSetsInt(int no) {
		set = new int [no];
	}

	public void clear() {
		int no = set.length;
		set = new int [no];
	}
	/**
	 * Union two disjoint sets using the height heuristic. root1 and root2 are
	 * distinct and represent set names.
	 * 
	 * @param root1
	 *            the root of set 1.
	 * @param root2
	 *            the root of set 2.
	 */
	public void union(int root1, int root2) {
		if (root1 == root2)
			return;

		int r1 = set[root1];
		int r2 = set[root2];

		if (r2 < r1) {
			// root2 is deeper
			// Make root2 new root
//			set.remove(root1);
//			set.put(root1, root2);
			set[root1] = root2;
		} else {
			if (r1 == r2) {
				// Update height if same
				r1--;
//				set.remove(root1);
//				set.put(root1, r1);
				set[root1] = r1;
			}
			// Make root1 new root
//			set.remove(root2);
//			set.put(root2, root1);
			set[root2] = root1;
		}
	}

	/**
	 * Perform a find with path compression.
	 * 
	 * @param x
	 *            the element being searched for.
	 * @return the set containing x.
	 */
	public int find(int x) {
		int f = set[x];
		if (f < 1) {
			return x;
		} else {
			int xx = find(f);
			set[x] = xx;
			return xx;
		}
	}


	// Test main; all finds on same output line should be identical
	public static void main(String[] args) {
		int numElements = 128;
		int numInSameSet = 16;

		// DisjointSets ds = new DisjointSets( numElements );
		DisjointSetsInt ds = new DisjointSetsInt(numElements + 1);
		int set1, set2;

		for (int k = 1; k < numInSameSet; k *= 2) {
			for (int j = 0; j + k < numElements; j += 2 * k) {
				set1 = ds.find(j + 1);
				set2 = ds.find(j + k + 1);
				ds.union(set1, set2);
			}
		}

		for (int i = 0; i < numElements; i++) {
			System.out.print(ds.find(i + 1) + "*");
			if (i % numInSameSet == numInSameSet - 1)
				System.out.println();
		}
		System.out.println();
	}
}
