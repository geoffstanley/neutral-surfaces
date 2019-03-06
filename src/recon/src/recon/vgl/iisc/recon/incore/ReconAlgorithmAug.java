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

import static vgl.iisc.recon.incore.TriangleDataPrim.getIndex;

import java.io.IOException;
import java.io.PrintStream;

import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.recon.incore.TriangleDataPrim.Vertex;
import vgl.iisc.utils.Utilities;

import vgl.iisc.reebgraph.ui.ReebGraphData; // GJS added

public class ReconAlgorithmAug {
	TriangleDataPrim data;
	StoreReebGraph rg;
	int noTris;
	
	int [] cpMap;
	DisjointSetsInt nodes;
	int [] nodeMap;
	ContourTree ct;
	MyArrays myArrays = new MyArrays();
	private boolean useAdj = true;

	public void run(double[] xx, double[] yy, double[] zz, double[] ff, int[] v1, int[] v2, int[] v3) { // GJS added
				
		useAdj = true; 		 // must be true for RAA.data to call addEdge, which is needed for making data.vertices[].link
		boolean storeCoords = false; // true results in the x,y,z vertex coordinates being stored in TriangleDataPrim java object
		
		int nV = xx.length; // Number of vertices
		int nF = v1.length; // Number of faces
		
		data = new TriangleDataPrim(useAdj, storeCoords, nV, nF);
		for( int i = 0; i < nV; i++ ) {
			data.addVertex(xx[i],yy[i],zz[i],ff[i]);
		}
		for( int i = 0; i < nF; i++ ) {
			data.addTriangle(v1[i],v2[i],v3[i]);
		}
		// Note: in MATLAB, should have subtracted 1 from faces so it references the vertex list with 0 indexing

		data.setupTriangles();
		setupData();
		orderVertices();
		computeReebGraph();

		rg.setup();
		rg.removeDeg2Nodes();
	}

	public class AugReeb { // GJS added
		public int nNodes;
		public int nArcs;
		public int[] node_v;
		public double[] node_fn;
		public byte[] node_type;
		public int[][] node_prev;
		public int[][] node_next;
		
		public int[] arc_from;
		public int[] arc_to;
		public int[] arc_id;

		public int[][] arc_segment;

		public AugReeb(int n, int a) {
			nNodes = n;
			nArcs = a;
			node_v = new int[n];
			node_fn = new double[n];
			node_type = new byte[n];
			node_prev = new int[n][];
			node_next = new int[n][];
			arc_from = new int[a];
			arc_to = new int[a];
			arc_id = new int[a];
		}
	}

	public AugReeb pack(int noVerts) { // GJS added
		return pack(noVerts, 0);
	}
	
	public AugReeb pack(int noVerts, int shift) { // GJS added
		
		ReebGraphData rgd = rg.rg.returnReebGraph();

		AugReeb ar = new AugReeb(rgd.nodes.length, rgd.arcs.length);

		// Iterate over arcs. 
		for( int i = 0; i < rgd.arcs.length; i++ ) {
			ar.arc_from[i] = rgd.arcs[i].from + shift;
			ar.arc_to[i] = rgd.arcs[i].to + shift;
			ar.arc_id[i] = rgd.arcs[i].id + shift;
		}

		ar.arc_segment = rg.rg.returnSegments(noVerts, shift);
		
		// Iterate over nodes:
		for( int i = 0; i < rgd.nodes.length; i++ ) {
			ar.node_v[i] = rgd.nodes[i].v + shift;
			ar.node_fn[i] = rgd.nodes[i].fn;
			ar.node_type[i] = rgd.nodes[i].type;
			
			ar.node_prev[i] = new int[rgd.nodes[i].prev.size()];
			for( int j = 0; j < ar.node_prev[i].length; j++ ) {
				ar.node_prev[i][j] = rgd.nodes[i].prev.get(j) + shift;
			}

			ar.node_next[i] = new int[rgd.nodes[i].next.size()];
			for( int j = 0; j < ar.node_next[i].length; j++ ) {
				ar.node_next[i][j] = rgd.nodes[i].next.get(j) + shift;
			}
			
		}
		
		return ar;

	}

	
	public void computeReebGraph(MeshLoader loader, String ftype) throws IOException {
		System.out.println("Loading mesh .... ");
		data = new TriangleDataPrim(useAdj);
		long loadTime = System.nanoTime();
		data.loadData(loader, ftype);
		loadTime = System.nanoTime() - loadTime;
		
		long gcTime = System.nanoTime();
		System.gc();System.gc();System.gc();System.gc();
		gcTime = System.nanoTime() - gcTime;
		
		System.out.println("Finished loading data .... ");
		double ct = System.nanoTime();
		setupData();
		orderVertices();
		computeReebGraph();
		long en = System.nanoTime();
		ct = (en - ct) / 1000000;
		
		System.out.println("No. of triangles : " + noTris);
		System.out.println("Time taken to compute Reeb Graph : " + ct + " ms");
	}
	
	double splitTime = 0;
	long splitTreeTime = 0;
	long joinTreeTime = 0;
	long mergeTime = 0;
	private void computeReebGraph() {
		nextV = data.noVertices;
		if((nextV & 1) == 1) {
			nextV ++;
		}

		findSplitTree();
		findJoinTree();
		ct.mergeTrees(rg);
	}

	void setupData() {
		maxStar = data.maxStar;
		noTris = data.triCt;
		
		ct = new ContourTree(data.noVertices, maxStar, noTris);
		
		lstarMap = new int[data.noVertices];
		ustarMap = new int[data.noVertices];
		criticalPts = new byte[data.noVertices];
		
		sv = new int[data.fnVertices.length];
		vMap = new int [data.fnVertices.length];
		starv = new int[data.noVertices];

		lstar = new Star(maxStar);
		ustar = new Star(maxStar);
		ustarTris = new int[maxStar];
		lstarTris = new int[maxStar];
		q = new int[noTris];
		/*
		 * For join and split trees 
		 */
		next = new int[maxStar];
		ntris = new int[maxStar];
		prev = new int[maxStar];
		oldVals = new int[maxStar];
		
		reached = new int[noTris];
		triangles = new int [noTris];
		jtriangles = new int [noTris];
		others = new boolean [noTris];
		

		for(int i = 0;i < noTris;i ++) {
			reached[i] = -1;
			triangles[i] = -1;
			jtriangles[i] = -1;
			if(i < data.fnVertices.length) {
				sv[i] = i;
				starv[i] = -1;
			}
		}
		
		// above code assumes that noTris >= noVertices. So handle that case
		
		if(noTris < data.fnVertices.length) {
			for(int i = noTris;i < data.fnVertices.length;i ++) {
				sv[i] = i;
				starv[i] = -1;
			}
		}
		link = new int[noTris];
		
		cpMap = new int[data.noVertices + noTris + noTris];
		nodes = new DisjointSetsInt(data.noVertices + noTris + noTris);
		nodeMap = new int[data.noVertices + noTris + noTris];
		
		
		nextSet = new int[data.noVertices + noTris + noTris];
		prevSet = new int[data.noVertices + noTris + noTris];
		
		for(int i = 0;i < nextSet.length;i ++) {
			nextSet[i] = -1;
			prevSet[i] = -1;
		}
		
		comps = new Component[maxStar];
		compNos = new int[maxStar];
		cps = new CriticalPoint[data.noVertices];
	}

	public void output(String op) {
		try {
			PrintStream p = new PrintStream(op);
			rg.setup();
			rg.removeDeg2Nodes();
			rg.outputReebGraph(p);
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void output(String op, String part) {
		try {
			PrintStream p = new PrintStream(op);
			PrintStream pt = new PrintStream(part);
			rg.setup();
			rg.removeDeg2Nodes();
			rg.outputReebGraph(p, pt, data.noVertices);
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	int [] sv;
	int [] vMap;
	int maxStar = 0;

	private void orderVertices() {
		myArrays.sort(sv);
		
		for(int i = 0;i < data.noVertices;i ++) {
			vMap[sv[i]] = i;
		}
	}

	public static final byte REGULAR = 0;
	public static final byte MINIMUM = 1;
	public static final byte SADDLE = 2;
	public static final byte MAXIMUM = 3;
	
	byte [] criticalPts;
	
	class Star {
		public Star(int maxStar) {
			adj = new int[2 * maxStar][2 * maxStar];
			ct = new int[2 * maxStar];
			tris = new int[2 * maxStar];
			reached = new boolean[2 * maxStar];
		}
		int [][] adj;
		int [] ct;
		boolean [] reached;
		int [] tris;
	}
	
	int [] starv;
	int [] lstarMap;
	int [] ustarMap;
	Star lstar;
	Star ustar;
	int [] ustarTris;
	int [] lstarTris;
	int ustarTrisct;
	int lstarTrisct;
	
	private byte getCriticality(int v) {
		int lct = 0;
		int uct = 0;
		for (int i = 0;i < data.vertices[v].star.length; i++) {
			int tt = data.vertices[v].star.array[i];
			int tin = tt * 3;
			int tv1 = data.triangles[tin];
			int tv2 = data.triangles[tin + 1];
			int tv3 = data.triangles[tin + 2];
			if (tv1 == v) {
				if(starv[tv2] != v) {
					starv[tv2] = v;
					ustarMap[tv2] = uct++;
					ustar.ct[ustarMap[tv2]] = 0;
					ustar.reached[ustarMap[tv2]] = false;
					ustar.tris[ustarMap[tv2]] = tt;
				}
				if(starv[tv3] != v) {
					starv[tv3] = v;
					ustarMap[tv3] = uct++;
					ustar.ct[ustarMap[tv3]] = 0;
					ustar.reached[ustarMap[tv3]] = false;
					ustar.tris[ustarMap[tv3]] = tt;
				}
				int v1 = ustarMap[tv2];
				int v2 = ustarMap[tv3];
				ustar.adj[v1][ustar.ct[v1]] = v2;
				ustar.adj[v2][ustar.ct[v2]] = v1;
				ustar.ct[v2] ++;
				ustar.ct[v1] ++;
			} else if (tv2 == v) {
				if(starv[tv1] != v) {
					starv[tv1] = v;
					lstarMap[tv1] = lct++;
					lstar.ct[lstarMap[tv1]] = 0;
					lstar.reached[lstarMap[tv1]] = false;
					lstar.tris[lstarMap[tv1]] = tt;
				}
				if(starv[tv3] != v) {
					starv[tv3] = v;
					ustarMap[tv3] = uct++;
					ustar.ct[ustarMap[tv3]] = 0;
					ustar.reached[ustarMap[tv3]] = false;
					ustar.tris[ustarMap[tv3]] = tt;
				}
			} else if (tv3 == v) {
				if(starv[tv1] != v) {
					starv[tv1] = v;
					lstarMap[tv1] = lct++;
					lstar.ct[lstarMap[tv1]] = 0;
					lstar.reached[lstarMap[tv1]] = false;
					lstar.tris[lstarMap[tv1]] = tt;
				}
				if(starv[tv2] != v) {
					starv[tv2] = v;
					lstarMap[tv2] = lct++;
					lstar.ct[lstarMap[tv2]] = 0;
					lstar.reached[lstarMap[tv2]] = false;
					lstar.tris[lstarMap[tv2]] = tt;
				}
				int v1 = lstarMap[tv1];
				int v2 = lstarMap[tv2];
				lstar.adj[v1][lstar.ct[v1]] = v2;
				lstar.adj[v2][lstar.ct[v2]] = v1;
				lstar.ct[v2] ++;
				lstar.ct[v1] ++;
			}
		}
		ustarTrisct = getStarComps(ustar, uct, ustarTris);
		lstarTrisct = getStarComps(lstar, lct, lstarTris);
		if(ustarTrisct == 0 && lstarTrisct == 0 || ustarTrisct == 1 && lstarTrisct == 1) {
			criticalPts[v] = REGULAR; 
			return REGULAR;
		}
		if(ustarTrisct == 0) {
			criticalPts[v] = MAXIMUM;
			return MAXIMUM;
		} 
		if(lstarTrisct == 0) {
			criticalPts[v] = MINIMUM;
			return MINIMUM;
		} else if(ustarTrisct > 1 || lstarTrisct > 1) {
			criticalPts[v] = SADDLE;
			return SADDLE;
		}
		return REGULAR;
	}
	
	int [] q;
	private int getStarComps(Star star, int ct, int [] starTris) {
		int comp = 0;
		for(int i = 0;i < ct;i ++) {
			if(!star.reached[i]) {
				starTris[comp ++] = star.tris[i];
				// do bfs
				int front = 0;
				int back = 0;
				q[back ++] = i;
				star.reached[i] = true;
				while(front < back) {
					int v = q[front ++];
					for(int j = 0;j < star.ct[v];j ++) {
						int vv = star.adj[v][j];
						if(!star.reached[vv]) {
							star.reached[vv] = true;
							q[back ++] = vv;
						}
					}
				}
			}
		}
		return comp;
	}
	
	
	/**
	 * 
	 * This is for computing join and split trees
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 */
	
	
	int [] reached;
	int [] triangles;
	int cpct = 1;
	int [] nextSet;
	int [] next;
	int [] oldVals;
	int nextct;
	
	int [] ntris;
	int ntrisct;
	
	int [] prevSet;
	int [] prev;
	int prevct;
	
	/* Begin Meta data*/
	int totCps = 0;
	int totLoops = 0;
	long unwantedSplit = 0;
	/* End Meta data*/
	
	int noCP = 0;
	int noLoops = 0;
	/* Split Tree */
	public void findSplitTree() {
		for(int i = data.noVertices - 1;i >= 0; i --) {
			int v = sv[i];
			byte cp = getCriticality(v);//criticalPts[v];

			if(cp == REGULAR) {
				processNormal(i, v);
			} else {
				totCps ++;
				noCP++;
				if(cp == MINIMUM) {
					if(ustarTrisct > 1) {
						if(processSaddle(i, v, true, true)) {
							noCP += cps[i].trisct;
							noLoops += cps[i].trisct;
							processAdditionalVertices(i);
							processSaddle(i, v, true, false);
						}
					} else {
						processSaddle(i, v, false, false);	
					}
				} else if(cp == MAXIMUM) {
					processMaximum(i, v);
				} else {
					if(ustarTrisct > 1) {
						if(processSaddle(i, v, true, true)) {
							noCP += cps[i].trisct;
							noLoops += cps[i].trisct;
							processAdditionalVertices(i);
							processSaddle(i, v, true, false);
						}
					} else if (lstarTrisct > 1) {
						processSaddle(i, v, false, false);
					} else {
						Utilities.er("Shouldn't come here!!!");
					}
				}
			}
		}
	}
	
	void processMinimum(int i,int v) {
		nextct = 0;
		Vertex vv;
		vv = data.vertices[v];
		for(int x = 0;x < vv.star.length; x++) {
			int tin = vv.star.array[x];
			if(triangles[tin] != -1) {
				int comp = nodes.find(triangles[tin]);
				if(nextSet[comp] != i) {
					nextSet[comp] = i;
					next[nextct ++] = comp;	
				}
			} else {
				Utilities.er("Should not come here");
			}
		}
		
		if(nextct > 0) {
			for(int x = 0;x < nextct;x ++) {
				int comp = nodes.find(next[x]);
				int p = cpMap[nodeMap[comp]];
				ct.addSplitArc(v, p);
			}
		}
	}
	
	void processMaximum(int i,int v) {
		Vertex vv;
		vv = data.vertices[v];

		cpMap[cpct] = v;
		for(int x = 0;x < vv.star.length; x++) {
			int tin = vv.star.array[x];
			triangles[tin] = cpct;
		}
		int comp = nodes.find(cpct);
		nodeMap[comp] = cpct;
		cpct ++;
	}

	void processNormal(int i,int v) {
		nextct = 0;
		ntrisct = 0;
		
		for(int x = 0;x < data.vertices[v].star.length; x++) {
			int tin = data.vertices[v].star.array[x];
			if(triangles[tin] != -1) {
				int comp = nodes.find(triangles[tin]);
				if(nextSet[comp] != i) {
					nextSet[comp] = i;
					next[nextct ++] = comp;	
				}
			}
			int tv3 = data.triangles[tin*3 + 2];
			if(tv3 == v) {
				ntris[ntrisct ++] = tin;
			}
		}
		if(nextct > 1) {
			Utilities.er("Regular vertex cannot have two nexts !!");
		}
		if(nextct == 1) {			
			cpMap[cpct] = v;
			int marker = cpct;
					
			for(int x = 0;x < ntrisct;x ++) {
				int tin = ntris[x];
				triangles[tin] = marker;
			}
			int comp = nodeMap[next[0]];
			int p = cpMap[comp];
			nodes.union(nodes.find(comp), nodes.find(marker));
			ct.addSplitArc(v, p);
			nodeMap[nodes.find(marker)] = cpct;
			cpct ++;
		}
	}
	int noSplits = 0;
	long splittimes = 0;
	
	Component [] comps;
	int [] compNos;
	
	CriticalPoint cp = null;
	
	boolean processSaddle(int i, int v, boolean split, boolean surgery) {
		nextct = 0;
		ntrisct = 0;
		
		Vertex vv;
		vv = data.vertices[v];
		
		for(int x = 0;x < vv.star.length; x++) {
			int tin = vv.star.array[x];
			if(triangles[tin] != -1) {
				int comp = nodes.find(triangles[tin]);
				if(nextSet[comp] != i) {
					oldVals[nextct] = nextSet[comp];
					nextSet[comp] = i;
					next[nextct ++] = comp;	
				}
			}
			int tv3 = data.triangles[tin * 3 + 2];
			if(tv3 == v) {
				ntris[ntrisct ++] = tin;
			}
		}
		if(nextct > 0) {
			if(split && surgery) {
				int pv = noActualSplits;
				cp = null;
				/* begin for meta data purpose*/
				long tm = 0; 
				/* end for meta data purpose*/
				if(nextct < ustarTrisct) {
					
					/* begin for meta data purpose*/
					totLoops ++;
					/* end for meta data purpose*/
					
					long st = System.nanoTime();
					int cct = 0;
					for(int x = 0;x < ustarTrisct; x ++) {
						int tr = ustarTris[x];
						int c = nodes.find(triangles[tr]);
						Component cc = null;
						for(int y = 0;y < cct;y ++) {
							if(compNos[y] == c) {
								cc = comps[y];
								break;
							}
						}
						if(cc == null) {
							cc = new Component();
							comps[cct] = cc;
							compNos[cct ++] = c;
						}
						cc.tris[cc.trisct ++] = (tr);
					}
					for (int x = 0;x < cct;x ++) {
						Component comp = comps[x];
						if (comp.trisct > 1) {
							split(i, comp.tris, comp.trisct);
							noSplits++;
						}
					}
					long en = System.nanoTime();
					
					splittimes += (en - st);
					tm = (en - st);
				}
				if(pv != noActualSplits) {
					// split occured
					cps[i] = cp;
					// revert currently changed component info.
					for(int x = 0;x < nextct;x ++) {
						int comp = next[x];
						nextSet[comp] = oldVals[x];
					}
					
					return true;
				}
				/* begin for meta data purpose*/
				else {
					unwantedSplit += tm;
				}
				/* end for meta data purpose*/
			}
			cpMap[cpct] = v;
			int marker = cpct;
					
			for(int x = 0;x < ntrisct;x ++) {
				int tin = ntris[x];
				triangles[tin] = marker;
			}
			for(int x = 0;x < nextct;x ++) {
				int comp = nodeMap[next[x]];
				int p = cpMap[comp];
				nodes.union(nodes.find(comp), nodes.find(marker));
				ct.addSplitArc(v, p);
			}
			nodeMap[nodes.find(marker)] = cpct;
			cpct ++;
		} else {
			Utilities.er("Shoudn't come here !!");
		}
		
		return false;
	}
	
	class Component {
		int[] tris = new int[maxStar];
		int trisct = 0;
	}
	int nextV;
	int noActualSplits = 0;
	boolean [] others;
	int oSize;
	int [] link;
	int linkct = 0;
	
	class CriticalPoint {
		int [] tris;
		int [] comp;
		int trisct = 0;
	}
	
	CriticalPoint [] cps;

	private void split(int i, int [] tris, int trisct) {
		oSize = trisct;
		for(int x = 0;x < trisct;x ++) {
			others[tris[x]] = true;
		}
		int ct = nextV;
		int first = ct;
		int no = 0;
		int tot = trisct;
		for(int xx = 0;xx < trisct;xx ++) {
			linkct = 0;
			int pv = nextV;
			nextV += 2;
			int compt = tris[xx];
			if(!others[compt]) {
				nextV -= 2;
				continue;
			}
			others[compt] = false;
			oSize --;
			int front = 0;
			int back = 0;
			
			q[back++] = compt;
			reached[compt] = ct;
			while(front < back) {
				int tin = q[front ++];
				if(others[tin]) {
					others[tin] = false;
					tot --;
					oSize --;
					if(oSize == 0 && ct == first) {
						// just one component
						return;
					}
				}
				
				int tIn = tin * 3;
				int tv1 = data.triangles[tIn];
				int tv2 = data.triangles[tIn + 1];
				int tv3 = data.triangles[tIn + 2];

				reached[tin] = ct;
				link[linkct ++] = tin;
				if(vMap[tv1] <= i && vMap[tv3] > i) {
					if(vMap[tv2] > i) {
						// v1 v2 and v1 v3
						if(!useAdj) {
							MyIntList fan = data.getTriangleFan(tv1, tv2);
							for(int fi = 0;fi < fan.length;fi ++) {
								int tt = fan.array[fi];
								if(tt != -1 && reached[tt] != ct) {
									q[back ++] = tt;
									reached[tt] = ct;
								}
							}
						} else {
							int tt = data.adjTriangle[getIndex(tin,0,0)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
							tt = data.adjTriangle[getIndex(tin,0,1)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					} else {
						// v2 v3 and v1 v3
						if(!useAdj) {
							MyIntList fan = data.getTriangleFan(tv2, tv3);
							for(int fi = 0;fi < fan.length;fi ++) {
								int tt = fan.array[fi];
								if(tt != -1 && reached[tt] != ct) {
									q[back ++] = tt;
									reached[tt] = ct;
								}
							}
						} else {
							int tt = data.adjTriangle[getIndex(tin,1,0)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
							tt = data.adjTriangle[getIndex(tin,1,1)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					}
					if(!useAdj) {
						MyIntList fan = data.getTriangleFan(tv1, tv3);
						for(int fi = 0;fi < fan.length;fi ++) {
							int tt = fan.array[fi];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					} else {
						int tt = data.adjTriangle[getIndex(tin,2,0)];
						if(tt != -1 && reached[tt] != ct) {
							q[back ++] = tt;
							reached[tt] = ct;
						}
						tt = data.adjTriangle[getIndex(tin,2,1)];
						if(tt != -1 && reached[tt] != ct) {
							q[back ++] = tt;
							reached[tt] = ct;
						}
					}
				}
			}
			no ++;
			if(cp == null) {
				cp = new CriticalPoint();
				cp.tris = new int[ustarTrisct];
				cp.comp = new int[ustarTrisct];
			}
			cp.comp[cp.trisct] = pv; // even
			cp.tris[cp.trisct ++] = compt;
			ct = nextV;
			if(no == tot - 1) {
				// we have the required no. of plsits
				while(xx < trisct) {
					others[tris[xx++]] = false;
				}
				break;
			}
		}
		noActualSplits += no;
	}

	private void processAdditionalVertices(int i) {
		CriticalPoint cp = cps[i];
		// odd minimum, even maximum
		for(int x = cp.trisct - 1;x >= 0;x --) {
			// The two new vertices
			int max = cp.comp[x]; // even
			int min = max + 1;
			// first process the min, then the max
			
			// For min
			cpMap[cpct] = min;
			int marker = cpct;
			cpct ++;		
			
			// For max
			cpMap[cpct] = max;
			
			
			// other processing
			nextct = 0;
			ntrisct = 0;
			
			int front = 0;
			int back = 0;
			int compt = cp.tris[x];
			q[back++] = compt;
			int ct = min;
			reached[compt] = ct;
			
			
			while(front < back) {
				int tin = q[front ++];
				
				int tIn = tin * 3;
				int tv1 = data.triangles[tIn];
				int tv2 = data.triangles[tIn + 1];
				int tv3 = data.triangles[tIn + 2];

				reached[tin] = ct;
				// process triangle tin as it is to be
				// first min
				
				if(triangles[tin] != -1) {
					int comp = nodes.find(triangles[tin]);
					if(nextSet[comp] != ct) {
						nextSet[comp] = ct;
						next[nextct ++] = comp;	
					}
				}
				triangles[tin] = cpct;
				
				if(vMap[tv1] <= i && vMap[tv3] > i) {
					if(vMap[tv2] > i) {
						// v1 v2 and v1 v3
						if(!useAdj) {
							MyIntList fan = data.getTriangleFan(tv1, tv2);
							for(int fi = 0;fi < fan.length;fi ++) {
								int tt = fan.array[fi];
								if(tt != -1 && reached[tt] != ct) {
									q[back ++] = tt;
									reached[tt] = ct;
								}
							}
						} else {
							int tt = data.adjTriangle[getIndex(tin,0,0)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
							tt = data.adjTriangle[getIndex(tin,0,1)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					} else {
						// v2 v3 and v1 v3
						if(!useAdj) {
							MyIntList fan = data.getTriangleFan(tv2, tv3);
							for(int fi = 0;fi < fan.length;fi ++) {
								int tt = fan.array[fi];
								if(tt != -1 && reached[tt] != ct) {
									q[back ++] = tt;
									reached[tt] = ct;
								}
							}
						} else {
							int tt = data.adjTriangle[getIndex(tin,1,0)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
							tt = data.adjTriangle[getIndex(tin,1,1)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					}
					if(!useAdj) {
						MyIntList fan = data.getTriangleFan(tv1, tv3);
						for(int fi = 0;fi < fan.length;fi ++) {
							int tt = fan.array[fi];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					} else {
						int tt = data.adjTriangle[getIndex(tin,2,0)];
						if(tt != -1 && reached[tt] != ct) {
							q[back ++] = tt;
							reached[tt] = ct;
						}
						tt = data.adjTriangle[getIndex(tin,2,1)];
						if(tt != -1 && reached[tt] != ct) {
							q[back ++] = tt;
							reached[tt] = ct;
						}
					}
				}
			}

			for(int xx = 0;xx < nextct;xx ++) {
				int comp = nodeMap[next[xx]];
				int p = cpMap[comp];
				nodes.union(nodes.find(comp), nodes.find(marker));
				this.ct.addSplitArc(min, p);
			}
			nodeMap[nodes.find(marker)] = marker;
			
			// for max
			int comp = nodes.find(cpct);
			nodeMap[comp] = cpct;
			cpct ++;
		}
	}
	
	/* Join Tree */
	int [] jtriangles;
	public void findJoinTree() {
		rg = new StoreReebGraph(cpct);
		nodes.clear();
		cpct = 1;
		for(int i = 0;i < sv.length;i ++) {
			int v = sv[i];
			
			if(criticalPts[v] == REGULAR) {
				rg.addNode(v, data.fnVertices[v], SADDLE);
				processNormalJ(i, v);
			} else {
				byte cp = criticalPts[v];
				rg.addNode(v, data.fnVertices[v], cp);
				if(cp == SADDLE) {
					processSaddleJ(i, v);
					if(cps[i] != null) {
						processAdditionalVerticesJ(i);
					}
				} else if(cp == MINIMUM) {
					processMinimumJ(i, v);
					if(cps[i] != null) {
						processAdditionalVerticesJ(i);
					}
				} else {
					processSaddleJ(i, v);
				}
			}
		}
	}
	
	void processMinimumJ(int i, int v) {
		cpMap[cpct] = v;
		Vertex vv;
		vv = data.vertices[v];

		for(int x = 0;x < vv.star.length; x++) {
			int tin = vv.star.array[x];
			jtriangles[tin] = cpct;
		}
		int comp = nodes.find(cpct);
		nodeMap[comp] = cpct;
		cpct ++;
	}
	
	void processMaximumJ(int i, int v) {
		prevct = 0;
		Vertex vv;
		vv = data.vertices[v];

		for(int x = 0;x < vv.star.length; x++) {
			int tin = vv.star.array[x];
			if(jtriangles[tin] != -1) {
				int comp = nodes.find(jtriangles[tin]);
				if(prevSet[comp] != i) {
					prevSet[comp] = i;
					prev[prevct ++] = comp;
				}
			} else {
				Utilities.er("Shoul not come here");
			}
		}
		
		if(prevct > 0) {
			for(int x = 0;x < prevct;x ++) {
				int comp = nodes.find(prev[x]);
				int p = cpMap[nodeMap[comp]];
				ct.addJoinArc(p, v);
			}
		}
	}

	void processNormalJ(int i, int v) {
		prevct = 0;
		ntrisct = 0;
		
		for(int x = 0;x < data.vertices[v].star.length; x++) {
			int tin = data.vertices[v].star.array[x];
			if(jtriangles[tin] != -1) {
				int comp = nodes.find(jtriangles[tin]);
				if(prevSet[comp] != i) {
					prevSet[comp] = i;
					prev[prevct ++] = comp;
				}
			}
			int tv1 = data.triangles[tin * 3];
			if(tv1 == v) {
				ntris[ntrisct ++] = tin;
			}
		}
		if(prevct > 1) {
			Utilities.er("Regular vertex cannot have two prevs !!");
		}
		if(prevct == 1) {			
			cpMap[cpct] = v;
			int marker = cpct;
					
			for(int x = 0;x < ntrisct;x ++) {
				int tin = ntris[x];
				jtriangles[tin] = marker;
			}
			
			if(prevct > 0) {
				for(int x = 0;x < prevct;x ++) {
					int comp = nodeMap[prev[x]];
					int p = cpMap[comp];
					nodes.union(nodes.find(comp), nodes.find(marker));
					ct.addJoinArc(p, v);
				}
			}
			nodeMap[nodes.find(marker)] = cpct;
			cpct ++;
		}
	}
	
	void processSaddleJ(int i, int v) {
		prevct = 0;
		ntrisct = 0;
		Vertex vv;
		vv = data.vertices[v];
		for(int x = 0;x < vv.star.length; x++) {
			int tin = vv.star.array[x];
			if(jtriangles[tin] != -1) {
				int comp = nodes.find(jtriangles[tin]);
				if(prevSet[comp] != i) {
					prevSet[comp] = i;
					prev[prevct ++] = comp;
				}
			}
			int tv1 = data.triangles[tin * 3];
			if(tv1 == v) {
				ntris[ntrisct ++] = tin;
			}
		}
		cpMap[cpct] = v;
		int marker = cpct;
				
		for(int x = 0;x < ntrisct;x ++) {
			int tin = ntris[x];
			jtriangles[tin] = marker;
		}
		
		if(prevct > 0) {
			for(int x = 0;x < prevct;x ++) {
				int comp = nodeMap[prev[x]];
				int p = cpMap[comp];
				nodes.union(nodes.find(comp), nodes.find(marker));
				ct.addJoinArc(p, v);
			}
		}
		nodeMap[nodes.find(marker)] = cpct;
		cpct ++;
	}

	private void processAdditionalVerticesJ(int i) {
		CriticalPoint cp = cps[i];
		// odd minimum, even maximum
		for(int x = 0;x < cp.trisct;x ++) {
			// The two new vertices
			int max = cp.comp[x]; // even
			int min = max + 1;
			
			rg.addNode(max, 0, (byte)0);
//			rg.addNode(min, 0, (byte)0);
			
			// first process the max, then the min
			// For max
			cpMap[cpct] = max;
			int marker = cpct;
			cpct ++;		
			
			// For min
			cpMap[cpct] = min;
			
			
			// other processing
			prevct = 0;
			ntrisct = 0;
			
			int front = 0;
			int back = 0;
			int compt = cp.tris[x];
			q[back++] = compt;
			int ct = max;
			reached[compt] = ct;
			
			
			while(front < back) {
				int tin = q[front ++];
				
				int tIn = tin * 3;
				int tv1 = data.triangles[tIn];
				int tv2 = data.triangles[tIn + 1];
				int tv3 = data.triangles[tIn + 2];

				reached[tin] = ct;
				// process triangle tin as it is to be
				// first max
				
				if(jtriangles[tin] != -1) {
					int comp = nodes.find(jtriangles[tin]);
					if(prevSet[comp] != ct) {
						prevSet[comp] = ct;
						prev[prevct ++] = comp;	
					}
				}
				jtriangles[tin] = cpct;
				
				if(vMap[tv1] <= i && vMap[tv3] > i) {
					if(vMap[tv2] > i) {
						// v1 v2 and v1 v3
						if(!useAdj) {
							MyIntList fan = data.getTriangleFan(tv1, tv2);
							for(int fi = 0;fi < fan.length;fi ++) {
								int tt = fan.array[fi];
								if(tt != -1 && reached[tt] != ct) {
									q[back ++] = tt;
									reached[tt] = ct;
								}
							}
						} else {
							int tt = data.adjTriangle[getIndex(tin,0,0)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
							tt = data.adjTriangle[getIndex(tin,0,1)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					} else {
						// v2 v3 and v1 v3
						if(!useAdj) {
							MyIntList fan = data.getTriangleFan(tv2, tv3);
							for(int fi = 0;fi < fan.length;fi ++) {
								int tt = fan.array[fi];
								if(tt != -1 && reached[tt] != ct) {
									q[back ++] = tt;
									reached[tt] = ct;
								}
							}
						} else {
							int tt = data.adjTriangle[getIndex(tin,1,0)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
							tt = data.adjTriangle[getIndex(tin,1,1)];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					}
					if(!useAdj) {
						MyIntList fan = data.getTriangleFan(tv1, tv3);
						for(int fi = 0;fi < fan.length;fi ++) {
							int tt = fan.array[fi];
							if(tt != -1 && reached[tt] != ct) {
								q[back ++] = tt;
								reached[tt] = ct;
							}
						}
					} else {
						int tt = data.adjTriangle[getIndex(tin,2,0)];
						if(tt != -1 && reached[tt] != ct) {
							q[back ++] = tt;
							reached[tt] = ct;
						}
						tt = data.adjTriangle[getIndex(tin,2,1)];
						if(tt != -1 && reached[tt] != ct) {
							q[back ++] = tt;
							reached[tt] = ct;
						}
					}
				}
			}

			for(int xx = 0;xx < prevct;xx ++) {
				int comp = nodeMap[prev[xx]];
				int p = cpMap[comp];
				nodes.union(nodes.find(comp), nodes.find(marker));
				this.ct.addJoinArc(p, max);
			}
			nodeMap[nodes.find(marker)] = marker;
			
			// for min
			int comp = nodes.find(cpct);
			nodeMap[comp] = cpct;
			cpct ++;
		}
	}
	
	public class MyArrays {

		private static final int INSERTIONSORT_THRESHOLD = 7;

		public void sort(int [] a) {
			int [] aux = clone(a);
			mergeSort(aux, a, 0, a.length, 0);
		}
		
		private int [] clone(int [] a) {
			int[] aux = new int[a.length];
			for(int i = 0;i < a.length;i ++) {
				aux[i] = a[i];
			}
			return aux;
		}
		private void mergeSort(int[] src, int[] dest, int low, int high, int off) {
			int length = high - low;

			// Insertion sort on smallest arrays
			if (length < INSERTIONSORT_THRESHOLD) {
				for (int i = low; i < high; i++)
					for (int j = i; j > low && compare(dest[j - 1], dest[j]) > 0; j--)
						swap(dest, j, j - 1);
				return;
			}

			// Recursively sort halves of dest into src
			int destLow = low;
			int destHigh = high;
			low += off;
			high += off;
			int mid = (low + high) >>> 1;
			mergeSort(dest, src, low, mid, -off);
			mergeSort(dest, src, mid, high, -off);

			// If list is already sorted, just copy from src to dest. This is an
			// optimization that results in faster sorts for nearly ordered lists.
			if (compare(src[mid - 1], src[mid]) <= 0) {
				System.arraycopy(src, low, dest, destLow, length);
				return;
			}

			// Merge sorted halves (now in src) into dest
			for (int i = destLow, p = low, q = mid; i < destHigh; i++) {
				if (q >= high || p < mid && compare(src[p], src[q]) <= 0)
					dest[i] = src[p++];
				else
					dest[i] = src[q++];
			}
		}

		private void swap(int[] x, int a, int b) {
			int t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}

	public int compare(int o1, int o2) {
		if(data.fnVertices[o1] < data.fnVertices[o2] || (data.fnVertices[o1] == data.fnVertices[o2] && o1 < o2)) {
			return -1;
		}
		return 1;
	}

	public void useAdjacencies(boolean adj) {
		useAdj = adj;
	}

}
