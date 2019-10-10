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

import static vgl.iisc.utils.Utilities.er;

import java.util.ArrayList;
import java.util.Arrays;

import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.external.types.Simplex;

public class TriangleDataPrim {

	public static final int PREV = 0;
	public static final int NEXT = 1;
	
	public class Vertex {
		public MyIntList star = new MyIntList();
		public MyIntList edges = new MyIntList();
	}

	public class Edge {
		int v1, v2;
	}

	public double[] fnVertices;
	public int noVertices;
	public int vertexCt;
	public Vertex[] vertices;

	MyIntList tris;
	public int [] triangles;
	MyIntList adjList;
	public int [] adjTriangle;
	public int triCt;

	public ArrayList<Edge> edges = new ArrayList<Edge>();
	public MyIntList nextTri = new MyIntList();
	public MyIntList firstTri = new MyIntList();
	int edCt;
	public int maxStar;
	
	boolean useAdj;
	boolean storeCoords;
	
	public double []x, y, z;
	
	public TriangleDataPrim(boolean useAdj) {
		this.useAdj = useAdj;
		this.storeCoords = false;
	}

	public TriangleDataPrim(boolean useAdj, boolean storeCoords) {
		this.useAdj = useAdj; 
		this.storeCoords = storeCoords;
	}

	public TriangleDataPrim(boolean useAdj, boolean storeCoords, int nv, int nt) { // GJS
		// For use with loading data from memory
		this.useAdj = useAdj; 
		this.storeCoords = storeCoords;
		setNoOfVertices(nv);
		this.tris = new MyIntList(nt);
		if(useAdj) {
			this.adjList = new MyIntList(nt * 6);
		}
	}
	
	public void addVertex(double xx, double yy, double zz, double ff) {
		if(storeCoords) {
			x[vertexCt] = xx;
			y[vertexCt] = yy;
			z[vertexCt] = zz;
		}
		vertices[vertexCt] = new Vertex();
		fnVertices[vertexCt++] = ff;
	}

	public void addTriangle(int v1, int v2, int v3) {
		int[] v = new int[3];
		if (less(v1, v2)) {
			if (less(v1, v3)) {
				if (less(v2, v3)) {
					v[0] = v1;
					v[1] = v2;
					v[2] = v3;
				} else {
					v[0] = v1;
					v[1] = v3;
					v[2] = v2;
				}
			} else {
				v[0] = v3;
				v[1] = v1;
				v[2] = v2;
			}
		} else {
			if (less(v2, v3)) {
				if (less(v1, v3)) {
					v[0] = v2;
					v[1] = v1;
					v[2] = v3;
				} else {
					v[0] = v2;
					v[1] = v3;
					v[2] = v1;
				}
			} else {
				v[0] = v3;
				v[1] = v2;
				v[2] = v1;
			}
		}
		v1 = v[0];
		v2 = v[1];
		v3 = v[2];

		if(isAdded(v1,v2,v3)) {
			return;
		}
		tris.add(v1);
		tris.add(v2);
		tris.add(v3);
		
		if(useAdj) {
			adjList.add(-1);
			adjList.add(-1);
			adjList.add(-1);
			adjList.add(-1);
			adjList.add(-1);
			adjList.add(-1);
			
			int ein = getPrev(v1, v2);
			if (ein == -1) {
				addNewEdge(v1, v2);
			} else {
				int preTri = nextTri.get(ein);
				nextTri.set(ein, triCt);
				addNextAdj(preTri, v1, v2);
				adjList.set(getIndex(triCt, 0,PREV), preTri);
			}
	
			ein = getPrev(v2, v3);
			if (ein == -1) {
				addNewEdge(v2, v3);
			} else {
				int preTri = nextTri.get(ein);
				nextTri.set(ein, triCt);
				addNextAdj(preTri, v2, v3);
				adjList.set(getIndex(triCt, 1,PREV), preTri);
			}
			
			ein = getPrev(v1, v3);
			if (ein == -1) {
				addNewEdge(v1, v3);
			} else {
				int preTri = nextTri.get(ein);
				nextTri.set(ein, triCt);
				addNextAdj(preTri, v1, v3);
				adjList.set(getIndex(triCt, 2,PREV), preTri);
			}
		}
		
		vertices[v1].star.add(triCt);
		vertices[v2].star.add(triCt);
		vertices[v3].star.add(triCt);
		maxStar = Math.max(maxStar, vertices[v1].star.length);
		maxStar = Math.max(maxStar, vertices[v2].star.length);
		maxStar = Math.max(maxStar, vertices[v3].star.length);
		triCt++;
	}

	
	private boolean isAdded(int v1, int v2, int v3) {
		int s1 = vertices[v1].star.size();
		int s2 = vertices[v2].star.size();
		int s3 = vertices[v3].star.size();

		if (s1 < s2 && s1 < s3) {
			return isPresent(v1, v1, v2, v3);
		} else if (s2 < s3) {
			return isPresent(v2, v1, v2, v3);
		} else {
			return isPresent(v3, v1, v2, v3);
		}
	}

	private boolean isPresent(int v, int v1, int v2, int v3) {
		int s = vertices[v].star.size();
		for (int i = 0; i < s; i++) {
			int t = vertices[v].star.get(i);
			t *= 3;
			int tv1 = tris.get(t);
			int tv2 = tris.get(t + 1);
			int tv3 = tris.get(t + 2);
			if (tv1 == v1 && tv2 == v2 && tv3 == v3) {
				return true;
			}
		}
		return false;
	}

	void addNextAdj(int tin, int v1, int v2) {
		int t = tin * 3;
		int tv1 = tris.get(t);
		int tv2 = tris.get(t + 1);
		
		if (tv1 == v1) {
			if (tv2 == v2) {
				// v1 v2
				adjList.set(getIndex(triCt, 0,NEXT), triCt);
			} else {
				// v1 v3
				adjList.set(getIndex(triCt, 2,NEXT), triCt);
			}
		} else {
			// v2 v3
			adjList.set(getIndex(triCt, 1,NEXT), triCt);
		}
	}

	private int getPrev(int v1, int v2) {
		int s1 = vertices[v1].edges.size();
		int s2 = vertices[v2].edges.size();
		int added = -1;
		if (s1 < s2) {
			added = isEdgePresent(v1, v1, v2);
		} else {
			added = isEdgePresent(v2, v1, v2);
		}
		return added;
	}

	private int isEdgePresent(int v, int v1, int v2) {
		int s = vertices[v].edges.size();

		for (int i = 0; i < s; i++) {
			int ein = vertices[v].edges.get(i);
			Edge ee = edges.get(ein);
			if (ee.v1 == v1 && ee.v2 == v2) {
				return ein;
			}
		}
		return -1;
	}

	void addNewEdge(int v1, int v2) {
		Edge e = new Edge();
		e.v1 = v1;
		e.v2 = v2;
		edges.add(e);
		nextTri.add(triCt);
		firstTri.add(triCt);
		vertices[v1].edges.add(edCt);
		vertices[v2].edges.add(edCt);
		edCt++;
	}

	private boolean less(int v1, int v2) {
		return ((fnVertices[v1] < fnVertices[v2]) || (fnVertices[v1] == fnVertices[v2] && v1 < v2));
	}

	public void loadData(MeshLoader loader, String ftype) {
		try {
			int nv = loader.getVertexCount();
			setNoOfVertices(nv);
			int nt = loader.getSimplexCount();
			tris = new MyIntList(nt);
			if(useAdj) {
				adjList = new MyIntList(nt * 6);
			}
			Simplex sim = loader.getNextSimplex();
			int f = Integer.parseInt(ftype);
			while (sim != null) {
				if (sim instanceof vgl.iisc.external.types.Vertex) {
					vgl.iisc.external.types.Vertex v = (vgl.iisc.external.types.Vertex) sim;
					double fn;
					if (f == 0) {
						fn = v.f;
					} else {
						fn = v.c[f - 1];
					}
					if(v.c.length < 3) {
						addVertex(v.c[0], v.c[1], 0, fn);	
					} else {
						addVertex(v.c[0], v.c[1], v.c[2], fn);
					}
				} else if (sim instanceof vgl.iisc.external.types.Triangle) {
					vgl.iisc.external.types.Triangle t = (vgl.iisc.external.types.Triangle) sim;
					addTriangle(t.v1, t.v2, t.v3);
				} else if (sim instanceof vgl.iisc.external.types.Tetrahedron) {
					vgl.iisc.external.types.Tetrahedron t = (vgl.iisc.external.types.Tetrahedron) sim;
					addTriangle(t.v1, t.v2, t.v3);
					addTriangle(t.v1, t.v2, t.v4);
					addTriangle(t.v1, t.v3, t.v4);
					addTriangle(t.v2, t.v3, t.v4);
				} else {
					er("Invalid Simplex");
				}
				sim = loader.getNextSimplex();
			}
//			pr("Finished reading data from file. Loading it......");
			setupTriangles();
//			pr("Successfully loaded Data");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	public void setupTriangles() {
		triangles = Arrays.copyOf(tris.array, tris.length);
		tris.clear();
		if(!useAdj) {
			return;
		}
		for(int i = 0;i < edCt;i ++) {
			int tin = nextTri.get(i);
			int fin = firstTri.get(i);
			Edge ein = edges.get(i);
			int ef = 0;
			int et = 0;
			int f = fin * 3;
			int e = tin * 3;
			
			int ev1 = triangles[e];
			int ev2 = triangles[e + 1];

			int fv1 = triangles[f];
			int fv2 = triangles[f + 1];
			
			if(ein.v1 == fv1) {
				if(ein.v2 == fv2) {
					ef = 0;
				} else {
					ef = 2;
				}
			} else {
				ef = 1;
			}

			if(ein.v1 == ev1) {
				if(ein.v2 == ev2) {
					et = 0;
				} else {
					et = 2;
				}
			} else {
				et = 1;
			}

			adjList.set(getIndex(fin, ef, PREV), tin);
			adjList.set(getIndex(tin, et, NEXT), fin);
		}
		adjTriangle = Arrays.copyOf(adjList.array, adjList.length);
		adjList.clear();
	}

	public void setNoOfVertices(int nv) {
		noVertices = nv;
		fnVertices = new double[noVertices];
		vertices = new Vertex[noVertices];
		
		if(storeCoords) {
			x = new double[nv];
			y = new double[nv];
			z = new double[nv];
		}
		vertexCt = 0;
		triCt = 0;
		edCt = 0;
	}
	
	MyIntList fan = new MyIntList();
	public MyIntList getTriangleFan(int v1, int v2) {
		fan.length = 0;
		int s1 = vertices[v1].star.length;
		int s2 = vertices[v2].star.length;
		if(s1 < s2) {
			return getFan(vertices[v1], v1, v2);
		} else {
			return getFan(vertices[v2], v1, v2);
		}
	}

	private MyIntList getFan(Vertex v, int v1, int v2) {
		for(int i = 0;i < v.star.length;i ++) {
			int t = v.star.array[i];
			int tin = t * 3;
			int tv1 = triangles[tin];
			int tv2 = triangles[tin + 1];
			int tv3 = triangles[tin + 2];
			if(tv1 == v1) {
				if(tv2 == v2 || tv3 == v2) {
					fan.add(t);
				}
			} else if(tv2 == v1) {
				if(tv3 == v2) {
					fan.add(t);
				}
			}
		}
		return fan;
	}
	
	public static int getIndex(int tin, int e, int n) {
		return (tin * 6 + e * 2 + n);
	}
}
