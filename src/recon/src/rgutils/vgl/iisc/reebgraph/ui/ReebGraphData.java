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
package vgl.iisc.reebgraph.ui;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

import static vgl.iisc.utils.Utilities.splitString;

public class ReebGraphData {
	public static final byte REGULAR = 0;
	public static final byte MINIMUM = 1;
	public static final byte SADDLE = 2;
	public static final byte MAXIMUM = 3;
	
	public static final String MINIMUM_STRING = "MINIMA";
	public static final String SADDLE_STRING = "SADDLE";
	public static final String MAXIMUM_STRING = "MAXIMA";
	
	public class Node {
		public int v;
		public double fn;
		public byte type;
		
		public ArrayList<Integer> prev = new ArrayList<Integer>();
		public ArrayList<Integer> next = new ArrayList<Integer>();
	}
	
	public class Arc {
		public int from;
		public int to;
		public int id;
		
		@Override
		public String toString() {
			String s = from + " " + to;
			return s;
		}
	}
	
	public int noNodes;
	public int noArcs;
	
	public Node [] nodes;
	public Arc [] arcs;
	public HashMap<Integer, Integer> nodeMap = new HashMap<Integer, Integer>();
	
	public ReebGraphData() {
		
	}
	
	public ReebGraphData(String file) {
		readReebGraph(file);
	}
	
	public ReebGraphData(int nv, int na) { // GJS
		noNodes = nv;
		noArcs = na;
		nodes = new Node[noNodes];
		arcs = new Arc[noArcs];
	}

	int curNode = 0; // GJS
	public void addNode(int v, double fn, String type) { // GJS
		nodes[curNode] = new Node();
		nodes[curNode].v = v;
		nodeMap.put(nodes[curNode].v, curNode);
		nodes[curNode].fn = fn;
		nodes[curNode].type = getType(type);
		curNode ++;
	}

	int curArc = 0; // GJS
	public void addArc(int v1, int v2) { // GJS
		arcs[curArc] = new Arc();
		arcs[curArc].id = curArc;
		arcs[curArc].from = nodeMap.get(v1);
		arcs[curArc].to = nodeMap.get(v2);
		nodes[arcs[curArc].from].next.add(curArc);
		nodes[arcs[curArc].to].prev.add(curArc);
		curArc ++;
	}

	public void readReebGraph(String file) {
		try {
			BufferedReader f = new BufferedReader(new FileReader(file));
			String s = f.readLine();
			String [] r = splitString(s);
			noNodes = Integer.parseInt(r[0].trim());
			noArcs = Integer.parseInt(r[1].trim());
			
			nodes = new Node[noNodes];
			arcs = new Arc[noArcs];
			
			
			for(int i = 0;i < noNodes;i ++) {
				s = f.readLine();
				r = splitString(s);
				nodes[i] = new Node();
				nodes[i].v = Integer.parseInt(r[0].trim());
				nodeMap.put(nodes[i].v, i);
				nodes[i].fn = Double.parseDouble(r[1].trim());
				nodes[i].type = getType(r[2].trim());
			}
			
			for(int i = 0;i < noArcs;i ++) {
				s = f.readLine();
				r = splitString(s);
				arcs[i] = new Arc();
				int v1 = Integer.parseInt(r[0].trim());
				int v2 = Integer.parseInt(r[1].trim());
				arcs[i].id = i;
				arcs[i].from = nodeMap.get(v1);
				arcs[i].to = nodeMap.get(v2);
				nodes[arcs[i].from].next.add(i);
				nodes[arcs[i].to].prev.add(i);
			}
			f.close();
		} catch(Exception e) {
			System.err.println("Error : Invalid file format");
			System.err.println(e.getMessage());
		}
	}

	public byte getType(String type) {
		if(MINIMUM_STRING.equalsIgnoreCase(type)) {
			return MINIMUM;
		}
		if(MAXIMUM_STRING.equalsIgnoreCase(type)) {
			return MAXIMUM;
		}
		if(SADDLE_STRING.equalsIgnoreCase(type)) {
			return SADDLE;
		}
		return REGULAR;
	}
	
	public String getTypeString(byte type) {
		if(MINIMUM == type) {
			return MINIMUM_STRING;
		}
		if(MAXIMUM == type) {
			return MAXIMUM_STRING;
		}
		if(SADDLE == type) {
			return SADDLE_STRING;
		}
		return "None";
	}
}
