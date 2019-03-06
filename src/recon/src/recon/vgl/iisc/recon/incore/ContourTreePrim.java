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

public class ContourTreePrim {
	
	class JSTree {
		int [] next;
		int [] prev;
		int [] nsize;
		int [] psize;
		boolean split;
		boolean [] present;
		public JSTree(int no, boolean split) {
			this.split = split;
			if(split) {
				next = new int[no * maxStar];
				nsize = new int[no];
				prev = new int[no];
				psize = new int[no];
			} else {
				prev = new int[no * maxStar];
				next = new int[no];
				psize = new int[no];
				nsize = new int[no];
			}
			present = new boolean[no];
		}
	}
	
	JSTree joinTree;
	JSTree splitTree;
	
	int nv; 
	int maxStar;
	int [] joinNodes;
	int [] splitNodes;
	int jct = 0;
	int sct = 0;
	
	public ContourTreePrim(int noVertices, int maxStar, int nt) {
		nv = noVertices;
		this.maxStar = maxStar / 2;

		joinTree = new JSTree(nv + nt, false);
		splitTree = new JSTree(nv + nt, true);
		
		joinNodes = new int[nv + nt];
		splitNodes = new int[nv + nt];
		
		q = new int[nv + nt];
	}

	public void addJoinArc(int from, int to) {
		if(!joinTree.present[from]) {
			joinTree.present[from] = true;
			joinNodes[jct ++] = from;
		}
		if(!joinTree.present[to]) {
			joinTree.present[to] = true;
			joinNodes[jct ++] = to;
		}
		
		joinTree.next[from] = to;
		joinTree.nsize[from] ++;
		joinTree.prev[to * maxStar + joinTree.psize[to]] = from;
		joinTree.psize[to] ++;
	}
	
	public void addSplitArc(int from, int to) {
		if(!splitTree.present[from]) {
			splitTree.present[from] = true;
			splitNodes[sct ++] = from;
		}
		if(!splitTree.present[to]) {
			splitTree.present[to] = true;
			splitNodes[sct ++] = to;
		}
		
		splitTree.next[from * maxStar + splitTree.nsize[from]] = to;
		splitTree.prev[to] = from;
		splitTree.nsize[from] ++;
		splitTree.psize[to] ++;
	}
	
	int [] q;
	public void mergeTrees(StoreReebGraph rg) {
		int front = 0;
		int back = 0;
		for(int x = 0;x < jct; x++) {
			int v = joinNodes[x];
			if(splitTree.nsize[v] + joinTree.psize[v] == 1) {
				q[back ++] = v;
			}
		}
		
		while(back > front + 1) {
			int xi = q[front ++];
			
			if(splitTree.nsize[xi] == 0 && splitTree.psize[xi] == 0) {
				if(!(joinTree.nsize[xi] == 0 && joinTree.psize[xi] == 0)) {
					System.out.println("Shoudnt happen!!!");
					System.exit(0);
				}
				continue;
			}
			if(splitTree.nsize[xi] == 0) {
				if(splitTree.psize[xi] > 1) {
					System.out.println("Can this happen too???");
					System.exit(0);
				}
				int xj = splitTree.prev[xi];
				remove(xi, joinTree);
				remove(xi, splitTree);
				int fr = xj;
				int to = xi;
				if(fr > nv) {
					if((fr & 1) == 1) {
						fr --;
					}
				}
				if(to > nv) {
					if((to & 1) == 1) {
						to --;
					}
				}
				rg.addArc(fr, to);
				if(splitTree.nsize[xj] + joinTree.psize[xj] == 1) {
					q[back ++] = xj;
				}
			} else {
				if(joinTree.nsize[xi] > 1) {
					System.out.println("Can this happen too???");
					System.exit(0);
				}
				if(joinTree.nsize[xi] == 0) {
					System.out.println("Can this happen too again???");
					System.exit(0);
				}
				int xj = joinTree.next[xi];
				
				remove(xi, joinTree);
				remove(xi, splitTree);
				
				int fr = xi;
				int to = xj;
				if(fr > nv) {
					if((fr & 1) == 1) {
						fr --;
					}
				}
				if(to > nv) {
					if((to & 1) == 1) {
						to --;
					}
				}
				rg.addArc(fr, to);

				if(splitTree.nsize[xj] + joinTree.psize[xj] == 1) {
					q[back ++] = xj;
				}
			}
		}
	}
	
	private void remove(int xi, JSTree tree) {
		if(tree.psize[xi] == 1 && tree.nsize[xi] == 1) {
			int p = tree.prev[xi];
			int n = tree.next[xi];
			int pmul = 1;
			int nmul = 1;
			if(tree.split) {
				n = tree.next[xi * maxStar];
				nmul = maxStar;
			} else {
				p = tree.prev[xi * maxStar];
				pmul = maxStar;
			}
			tree.psize[xi] = 0;
			tree.nsize[xi]  = 0;
			
			removeAndAdd(tree.next, p*nmul, tree.nsize[p], xi, n);
			removeAndAdd(tree.prev, n*pmul, tree.psize[n], xi, p);
		} else if(tree.psize[xi] == 0 && tree.nsize[xi] == 1) {
			int n = tree.next[xi];
			int pmul = 1;
			if(tree.split) {
				n = tree.next[xi * maxStar];
			} else {
				pmul = maxStar;
			}
			tree.nsize[xi] = 0;
			remove(tree.prev, n * pmul, tree.psize[n], xi);
			tree.psize[n] --;
		} else if(tree.psize[xi] == 1 && tree.nsize[xi] == 0) {
			int p = tree.prev[xi];
			int nmul = 1;
			if(!tree.split) {
				p = tree.prev[xi * maxStar];
			} else {
				nmul = maxStar;
			}
			tree.psize[xi] = 0;
			remove(tree.next, p * nmul, tree.nsize[p], xi);
			tree.nsize[p] --;
		} else {
			System.out.println("Can this too happen??????");
			System.exit(0);
		}
	}

	private void remove(int[] arr, int start, int arrSize, int xi) {
		for(int i = start;i < start + arrSize;i ++) {
			if(arr[i] == xi) {
				if(i != start + arrSize - 1) {
					arr[i] = arr[start + arrSize - 1];
				}
				return;
			}
		}
		System.out.println("Shouldn't happen");
		System.exit(0);
	}

	private void removeAndAdd(int[] arr, int start, int arrSize, int rem, int add) {
		for(int i = start;i < start + arrSize;i ++) {
			if(arr[i] == rem) {
				arr[i] = add;
				return;
			}
		}
		System.out.println("Shouldn't happen");
		System.exit(0);
	}
}
