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
 *	Modified by : Geoffrey J. Stanley
 *	Date : 25/11/2018
 *	Changes  : marked by GJS
 */
package vgl.iisc.recon.incore;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import vgl.iisc.utils.Utilities;
import vgl.iisc.reebgraph.ui.ReebGraphData; // GJS


public class CleanReebGraph implements Serializable {

	private static final long serialVersionUID = 1L;

	public static enum VertexType {
		MINIMA, MAXIMA, SADDLE, NONE
	}
	
	public class Arc implements Comparable<Arc>, Serializable {
		private static final long serialVersionUID = 1L;

		int v1;
		int v2;
		double fn;
		int icol;
		boolean removed;
//		ArrayList <Integer> comps = new ArrayList<Integer>();
		
		
//		ArrayList<Vertex> path = new ArrayList<Vertex>();
		HashSet<Integer> segment = new HashSet<Integer>();
		
		public boolean equals(Object obj) {
			Arc a = (Arc) obj;
			return (a.icol == icol);
		}

		public int compareTo(Arc o) {
			double ff = fn - o.fn;
			if(ff < 0)
				return -1;
			else if(ff > 0)
				return 1;
			return 0;
		}
	}
	
	public class Node implements Serializable {
		private static final long serialVersionUID = 1L;

		int v;
		boolean removed;
		double fn;
		VertexType type;
		ArrayList<Arc> prev = new ArrayList<Arc>();
		ArrayList<Arc> next = new ArrayList<Arc>();
	}
	
	public Node [] nodes;
	public Arc [] arcs;

	ArrayList<Node> an = new ArrayList<Node>();
	ArrayList<Arc> ar = new ArrayList<Arc>();
	public HashMap<Integer, Integer> vmap = new HashMap<Integer, Integer>();
	public int [] amap; // GJS
	int ct = 0;
	short ect = 0;
	double min = Double.MAX_VALUE;
	double max = -Double.MAX_VALUE;
	double persistence;
	
	public void addNode(int v, double fn, byte type) {
		Node n = new Node();
		n.fn = fn;
		n.v = v;
		if(type == ReConAlgorithmPrim.MINIMUM) {
			n.type = VertexType.MINIMA;
		}
		if(type == ReConAlgorithmPrim.MAXIMUM) {
			n.type = VertexType.MAXIMA;
		}
		if(type == ReConAlgorithmPrim.SADDLE) {
			n.type = VertexType.SADDLE;
		}
		an.add(n);
		vmap.put(v, ct);
		ct ++;
		
		max = Math.max(max, fn);
		min = Math.min(min, fn);
	}
	
	public void setup() {
		nodes = an.toArray(new Node[0]);
		arcs = ar.toArray(new Arc[0]);
		persistence = max - min;
	}
	
	public void addArc(int v1, int v2) {
		Arc a = new Arc();
		a.v1 = vmap.get(v1);
		a.v2 = vmap.get(v2);
		a.fn = an.get(a.v2).fn - an.get(a.v1).fn;

//		a.path = new ArrayList<Vertex>();
		a.segment.add(v1);
		a.segment.add(v2);
		ect ++;
		a.icol = ect;
		ar.add(a);
		
		an.get(a.v1).next.add(a);
		an.get(a.v2).prev.add(a);
	}

	public void removeDeg2Nodes() {
		// remove degree 2 vertices
		for(int i = 0;i < nodes.length;i ++) {
			if(nodes[i].prev.size() == 0 && nodes[i].next.size() == 0) {
				nodes[i].removed = true;
			}
			if(!nodes[i].removed && nodes[i].next.size() == 1 && nodes[i].prev.size() == 1) {
				mergeNode(i);
			}
		}
	}
	
	private void mergeNode(int i) {
		Arc e1 = nodes[i].prev.get(0);
		Arc e2 = nodes[i].next.get(0);
		
		nodes[i].removed = true;
		e1.v2 = e2.v2;
		e2.removed = true;
//		if(e2.path.size() != 0) {
//			e2.path.remove(0);
//		}
//		e1.path.addAll(e2.path);
		e1.segment.addAll(e2.segment);
		if(e1.icol > e2.icol) {
			e1.icol = e2.icol;
		}
		e1.fn += e2.fn;
		nodes[e1.v2].prev.remove(e2);
		nodes[e1.v2].prev.add(e1);
	}
	
	public void outputSegmentation(PrintStream p, int noVerts) {
		int [] part = new int [nodes.length];
		Arrays.fill(part, -1);
		for(int i = 0;i < arcs.length;i ++) {
			if(!arcs[i].removed) {
				for(Iterator<Integer> it = arcs[i].segment.iterator();it.hasNext();) {
					int v = it.next();
					if(v < noVerts) {
						part[v] = arcs[i].icol;
					}
				}
			}
		}
		p.println(noVerts);
		for(int i = 0;i < noVerts;i ++) {
			if(part[i] == -1) {
				Utilities.er("Vertex not part of arc!!");
			}
			p.println(part[i]);
		}
	}
	
	public int[][] returnSegments(int noVerts, int shift) { // GJS 
		// Turn the arcs segment HashMap into sorted int[][].
		// Also, only include vertices whose index is < noVerts. Can use this to exclude
		// the additional vertices, added in ReconAlgorithmAug/processAdditionalVerticesJ(),
		// that are made for the splitting into different regions without holes.
		int iA = 0;
		Arc a;
		int iV;
		for(int i = 0; i < arcs.length; i++) {
			if(!arcs[i].removed) {
				iA++ ;
			}
		}
		int[][] seg = new int[iA][];
		iA = 0;
		for(int i = 0; i < arcs.length; i++) {
			if(!arcs[i].removed) {
				a = arcs[i];
				iV = 0;
				for(Iterator<Integer> it = a.segment.iterator();it.hasNext();) {
					int v = it.next();
					if (v < noVerts) {
						iV++;
					}
				}
				seg[iA] = new int[iV];
				iV = 0;
				for(Iterator<Integer> it = a.segment.iterator();it.hasNext();) {
					int v = it.next();
					if (v < noVerts) {
						seg[iA][iV++] = v + shift;
					}
				}
				Arrays.sort(seg[iA]);
				iA++;
			}
		}
		return seg;
	}

	public void outputReebGraph(PrintStream p) {
		int nv = 0;
		int ne = 0;
		for(int i = 0;i < nodes.length;i ++) {
			if(!nodes[i].removed) {
				nv ++;
			}
		}
		
		for(int i = 0;i < arcs.length;i ++) {
			if(!arcs[i].removed) {
				arcs[i].icol = ne;
				ne ++;
			}
		}

		p.println(nv + " " + ne);
		System.out.println(nv + " " + ne);
		for(int i = 0;i < nodes.length;i ++) {
			if(!nodes[i].removed) {
				p.println(nodes[i].v + " " + nodes[i].fn + " " + nodes[i].type);
			}
		}
		
		for(int i = 0;i < arcs.length;i ++) {
			if(!arcs[i].removed) {
				p.print(nodes[arcs[i].v1].v + " " + nodes[arcs[i].v2].v + " ");
//				for(int j = 0;j < arcs[i].comps.size();j ++) {
//					p.print(arcs[i].comps.get(j) + " ");
//				}
				p.println();
			}
		}
	}
	
	public ReebGraphData returnReebGraph() { // GJS
		// Basically as in outputReebGraph but different output.		
		int nv = 0;
		for(int i = 0;i < nodes.length;i ++) {
			if(!nodes[i].removed) {
				nv ++;
			}
		}

		int ne = 0;
		for(int i = 0;i < arcs.length;i ++) {
			if(!arcs[i].removed) {
				arcs[i].icol = ne;
				ne ++;
			}
		}

		ReebGraphData rgd = new ReebGraphData(nv, ne);
		
		for(int i = 0;i < nodes.length;i ++) {
			if(!nodes[i].removed) {
				rgd.addNode(nodes[i].v, nodes[i].fn, nodes[i].type.toString());
			}
		}

		amap = new int[ne];
		ne = 0;
		for(int i = 0;i < arcs.length;i ++) {
			if(!arcs[i].removed) {
				amap[ne] = i;
				rgd.addArc( nodes[arcs[i].v1].v, nodes[arcs[i].v2].v );
				ne ++;
			}
		}

		return rgd;
	}

}
