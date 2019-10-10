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
package vgl.iisc.reebgraph.ui;

import java.util.ArrayList;
import java.util.Iterator;

import vgl.iisc.reebgraph.ui.BranchDecomposition.Branch;

public class SimplifyReebGraph {
	
	public class Partition {
		public ArrayList<Integer> children = new ArrayList<Integer>();
	}
	
	private BranchDecomposition decomp;
	private ReebGraphData rgData;
	
	int noBranches;
	int noNodes;
	double maxPersistence;

	public Partition [] edgePartition;
	public int lastIncludedBranch;
	public boolean [] drawNode;
	
	public SimplifyReebGraph(ReebGraphData rgData, BranchDecomposition decomp) {
		this.rgData = rgData;
		this.decomp = decomp;
		
		setup();
	}
	
	private void setup() {
		noBranches = decomp.branches.size();
		noNodes = rgData.noNodes;
		maxPersistence = decomp.branches.get(0).fn;
		
		edgePartition = new Partition[rgData.noArcs];
		
		for(int i = 0;i < rgData.noArcs;i ++) {
			edgePartition[i] = new Partition();
			edgePartition[i].children.add(i);
		}
		
		drawNode = new boolean[noNodes];
		lastIncludedBranch = noBranches - 1;
		for(int i = 0;i < noNodes;i ++) {
			drawNode[i] = true;
		}
	}
	
	public void simplify(double sim) {
		lastIncludedBranch = noBranches - 1;
		double val = sim * maxPersistence;
		for(int i = 0;i < rgData.noArcs;i ++) {
			edgePartition[i].children.clear();
			edgePartition[i].children.add(i);
		}
		
		for(int i = 0;i < noNodes;i ++) {
			drawNode[i] = false;
		}
		
		for(int i = noBranches-1;i >= 0;i --) {
			Branch br = decomp.branches.get(i);
			if(br.fn < val) {
				removeBranch(i, br);
				lastIncludedBranch --;
			} else {
				break;
			}
		}
		
		for(int i = 0;i <= lastIncludedBranch;i ++) {
			Branch br = decomp.branches.get(i);
			drawNode[br.from] = true;
			drawNode[br.to] = true;
		}
	}

	public void simplify(int simBranches) {
		lastIncludedBranch = noBranches - 1;
		for(int i = 0;i < rgData.noArcs;i ++) {
			edgePartition[i].children.clear();
			edgePartition[i].children.add(i);
		}
		
		for(int i = 0;i < noNodes;i ++) {
			drawNode[i] = false;
		}
		
		for(int i = noBranches-1;i >= 0;i --) {
			Branch br = decomp.branches.get(i);
			if(br.fn == 0 || lastIncludedBranch >= simBranches) {
				removeBranch(i, br);
				lastIncludedBranch --;
			} else {
				break;
			}
		}
		
		for(int i = 0;i <= lastIncludedBranch;i ++) {
			Branch br = decomp.branches.get(i);
			drawNode[br.from] = true;
			drawNode[br.to] = true;
		}
	}
	
	private void removeBranch(int in, Branch br) {
		Branch par = br.parent;
		if(decomp.nodeBranch[br.from] == par.id) {
			// top
			for(Iterator<Integer> it = rgData.nodes[br.from].next.iterator();it.hasNext();) {
				int e = it.next();
				if(par.arcs.contains(e)) {
					for(Iterator<Integer> at = br.arcs.iterator();at.hasNext();) {
						int ee = at.next();
						edgePartition[e].children.addAll(edgePartition[ee].children);
						edgePartition[ee].children.clear();
					}
				}
			}
		} else {
			// bottom
			for(Iterator<Integer> it = rgData.nodes[br.to].prev.iterator();it.hasNext();) {
				int e = it.next();
				if(par.arcs.contains(e)) {
					for(Iterator<Integer> at = br.arcs.iterator();at.hasNext();) {
						int ee = at.next();
						edgePartition[e].children.addAll(edgePartition[ee].children);
						edgePartition[ee].children.clear();
					}
				}
			}
		}
	}
}
