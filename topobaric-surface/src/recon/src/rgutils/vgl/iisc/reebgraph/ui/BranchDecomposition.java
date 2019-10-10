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

import static vgl.iisc.utils.Utilities.er;
import static vgl.iisc.utils.Utilities.pr;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.PriorityQueue;

public class BranchDecomposition {
    
	public enum BranchDrawType {
        BOTTOM, TOP, MIDDLE, LOOP
    };
    
	public class Branch implements Comparable<Branch> {
		public int from;
		public int to;
		public double fn;
		
		public ArrayList<Integer> arcs = new ArrayList<Integer>();
		
		public ArrayList<Branch> children = new ArrayList<Branch>();
		public Branch parent;
		
		public ArrayList<Branch> connectingChildren = new ArrayList<Branch>();
		public Branch anotherParent;
		public boolean connectingLoop = false;
		public boolean loop = false;
		public int id;
		
		public int compareTo(Branch o) {
			if(fn > o.fn) {
				return -1; 
			} else if(fn < o.fn) {
				return 1;
			}
			return 0;
		}
		
		@Override
		public String toString() {
			return id + ": " + from + " " + to;
		}
	}
	
	public ArrayList<Branch> branches = new ArrayList<Branch>();
	public Branch root;
	public int [] nodeBranch;
	public int [] arcBranch;
	public BranchDrawType [] type;
	private ReebGraphData data;
	private boolean [] used;
	// store prev edge
	private int [] prev;
	private int [] q;
	private int [] reached;
	private ArrayList<Integer> list = new ArrayList<Integer>();
	private int cpct = 1;
	private int curBranch = 0;
	private boolean first = true;
	
	double [] edgeWts;
	boolean persistence; 
	
	public BranchDecomposition(ReebGraphData rgData, boolean compute) {
		persistence = true;
		data = rgData;
		if(compute) {
			compute();
		}
	}
	
	public BranchDecomposition(ReebGraphData rgData, double [] edgeWts) {
		persistence = true;
		data = rgData;
		if(edgeWts != null) {
			this.edgeWts = edgeWts;
			persistence = false;
		}
		compute();
	}
	
	private void compute() {
		System.out.println("Creating branches");
		createBranches();
		System.out.println("Organizing branches");
		organizeBranches();
	}
	private void createBranches() {
		used = new boolean[data.noArcs];
		prev = new int[data.noNodes];
		reached = new int[data.noNodes];
		nodeBranch = new int[data.noNodes];
		q = new int[data.noNodes];
		type = new BranchDrawType[data.noArcs];
		arcBranch = new int [data.noArcs]; 
		for(int i = 0;i < data.noNodes;i ++) {
			reached[i] = -1;
			nodeBranch[i] = -1;
		}
		int remaining = data.noArcs;
		int old;
		firstTime = true;
		do {
			old = remaining;
			remaining = processNodes();
		} while(old > remaining);
		if(remaining != 0) {
			er("Bad Behaviour!!");
		}
	}
	
	boolean firstTime;
	private int processNodes() {
		int no = Math.max(data.noNodes / 100,1);
		int ret = 0;

		for(int i = 0;i < data.noNodes;i ++) {
			while(potentialNode(i)) {
				int pair = getPair(i);
				if(pair == -1 || pair == i) {
					ret++;
					break;
				}
				Branch br = new Branch();
				br.from = i;
				br.to = pair;
				if(persistence) {
					br.fn = data.nodes[pair].fn - data.nodes[i].fn;
				} else {
					br.fn = wt * (data.nodes[pair].fn - data.nodes[i].fn);
				}
				branches.add(br);
				br.id = curBranch;
				updateEdges(br);
				curBranch ++;
			}
			if(i % no == 0) {
//				System.out.println(i / no);
			}
			first = false;
		}
		return ret;
	}
	
	private int getPair(int st) {
		int front = 0;
		int back = 0;
		q[back ++] = st;
		list.clear();
		int flag = cpct;
		cpct ++;
		while(front < back) {
			int v = q[front ++];
			int proceed = 0;
			for(Iterator<Integer> it = data.nodes[v].next.iterator();it.hasNext();) {
				int e = it.next();
				if(!used[e] && reached[data.arcs[e].to] != flag && (nodeBranch[v] == -1 || v == st)) {
					q[back ++] = data.arcs[e].to;
					reached[data.arcs[e].to] = flag;
					proceed ++;
					prev[data.arcs[e].to] = e;
				}
			}
			if(proceed == 0) {
				// possible extremum
				if(first || nodeBranch[v] != -1 || nodeBranch[st] != -1) {
					list.add(v);
				}
			}
		}
		return getPairedVertex(st);
	}
	
	double wt = 0;
	private int getPairedVertex(int st) {
		int max = -1;
		wt = 0;
		for(Iterator<Integer> it = list.iterator();it.hasNext();) {
			int v = it.next();
			if(persistence) {
				if(max < v) {
					max = v;
				}
			} else {
				double nwt = getBranchWt(st, v);
				if(nwt > wt) {
					wt = nwt;
					max = v;
				}
			}
			if(!first && nodeBranch[v] == nodeBranch[st]) {
				// give preference to parallel loops
				max = v;
				if(!persistence) {
					wt = getBranchWt(st, v);
				}
				break;
			}
		}
		return max;
	}
	
	private double getBranchWt(int st, int en) {
		double wt = 0;
		int v = en;
		while(v != st) {
			int e = prev[v];
			wt += edgeWts[e];
			v = data.arcs[e].from;
		}
		return wt;
	}

	private void updateEdges(Branch br) {
		int v = br.to;
		if(nodeBranch[v] == -1) {
			nodeBranch[v] = curBranch;
		}
		while(v != br.from) {
			int e = prev[v];
			used[e] = true;
			br.arcs.add(0,e);
			v = data.arcs[e].from;
			if(nodeBranch[v] == -1) {
				nodeBranch[v] = curBranch;
			}
			
			int from = data.arcs[e].from;
			int to = data.arcs[e].to;
			
			if(nodeBranch[from] == nodeBranch[to]) {
				if(nodeBranch[from] == curBranch) {
					type[e] = BranchDrawType.MIDDLE;
				} else {
					type[e] = BranchDrawType.LOOP;
				}
			} else if(nodeBranch[from] == curBranch) {
				type[e] = BranchDrawType.TOP;
			} else if(nodeBranch[to] == curBranch) {
				type[e] = BranchDrawType.BOTTOM;
			}
		}
	}
	
	private boolean potentialNode(int i) {
		for(Iterator<Integer> it = data.nodes[i].next.iterator();it.hasNext();) {
			int e = it.next();
			if(!used[e]) {
				return true;
			}
		}
		return false;
	}
	
	private void organizeBranches() {
		System.out.println("Organizing branches");
		boolean foundRoot = false;
		for(Iterator<Branch> it = branches.iterator();it.hasNext();) {
			Branch br = it.next();
			if(data.nodes[br.from].type == ReebGraphData.MINIMUM && data.nodes[br.to].type == ReebGraphData.MAXIMUM) {
				// Root Branch
				if(foundRoot) {
					er("More than one root!!!!");
				}
				foundRoot = true;
				root = br;
			} else if(nodeBranch[br.from] == br.id) {
				Branch parent = branches.get(nodeBranch[br.to]);
				parent.children.add(br);
				br.parent = parent;
			} else if(nodeBranch[br.to] == br.id) {
				Branch parent = branches.get(nodeBranch[br.from]);
				parent.children.add(br);
				br.parent = parent;
			} else {
				if(nodeBranch[br.to] == nodeBranch[br.from]) {
//					pr("There's a loop!!");
					Branch parent = branches.get(nodeBranch[br.from]);
					parent.children.add(br);
					br.parent = parent;
					br.loop = true;
				} else {
//					pr("There's a connecting loop!!");
					br.connectingLoop = true;
					Branch parent = branches.get(nodeBranch[br.from]);
					br.parent = parent;
					parent.connectingChildren.add(br);
					
					parent = branches.get(nodeBranch[br.to]);
					br.anotherParent = parent;
					parent.connectingChildren.add(br);
				}
			}
		}
		
		pr("Ordering branches");
		int no = branches.size();
		branches.clear();
		PriorityQueue<Branch> pq = new PriorityQueue<Branch>();
		pq.add(root);
		boolean [] processed = new boolean[no];
		while(pq.size() > 0) {
			Branch br = pq.poll();
			branches.add(br);
			processed[br.id] = true;
			for(Iterator<Branch> it = br.children.iterator();it.hasNext();) {
				Branch ch = it.next();
				pq.add(ch);
			}
			for(Iterator<Branch> it = br.connectingChildren.iterator();it.hasNext();) {
				Branch ch = it.next();
				Branch par = null;
				if(ch.parent == br) {
					par = ch.anotherParent;
				} else {
					par = ch.parent;
				}
				if(processed[par.id]) {
					pq.add(ch);
				}
			}
		}
		// reorder ids;
		int [] newMap = new int[no];
		int ct = 0;
		for(Iterator<Branch> it = branches.iterator();it.hasNext();) {
			Branch br = it.next();
			newMap[br.id] = ct;
			br.id = ct ++;
			for(Iterator<Integer> ait = br.arcs.iterator();ait.hasNext();) {
				int a = ait.next();
				arcBranch[a] = br.id;
			}
		}
		for(int i = 0;i < nodeBranch.length;i ++) {
			nodeBranch[i] = newMap[nodeBranch[i]];
		}
		
		if(no != branches.size()) {
			System.out.println(no + " " + branches.size());
			er("Some branches are not accounted for!!");
		}
	}

	public void printBranches() {
		int ct = 0;
		for(Iterator<Branch> it = branches.iterator();it.hasNext();) {
			Branch br = it.next();
			System.out.println(br.from + " " + br.to);
			ct += br.arcs.size();
		}
		System.out.println("No. of Branches : " + branches.size());
		System.out.println("No. of arcs : " + ct + " " + data.noArcs);
		System.out.println("Root Branch : " + root.from + " " + root.to);
	}
	
	public static void main(String [] args) {
		String rg = "Data/test.rg";
		ReebGraphData data = new ReebGraphData();
		System.out.println("Loading data");
		data.readReebGraph(rg);
		System.out.println("Getting BD");
		BranchDecomposition bd = new BranchDecomposition(data, true);
		bd.printBranches();
	}
}
