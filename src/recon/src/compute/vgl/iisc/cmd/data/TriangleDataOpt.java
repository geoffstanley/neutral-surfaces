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
package vgl.iisc.cmd.data;

import static vgl.iisc.utils.Utilities.er;
import static vgl.iisc.utils.Utilities.pr;

import java.util.ArrayList;
import java.util.Iterator;

import vgl.iisc.external.loader.MeshLoader;
import vgl.iisc.external.types.Simplex;

public class TriangleDataOpt {

	public static final int PREV = 0;
	public static final int NEXT = 1;

	public class Triangle {

		int v1, v2, v3;
		int[][] adjTriangle = new int[3][2];

		public Triangle() {
			v1 = -1;
			v2 = -1;
			v3 = -1;
			adjTriangle[0][0] = -1;
			adjTriangle[0][1] = -1;
			adjTriangle[1][0] = -1;
			adjTriangle[1][1] = -1;
			adjTriangle[2][0] = -1;
			adjTriangle[2][1] = -1;
		}
	}

	public class Vertex {
		ArrayList<Integer> starList = new ArrayList<Integer>();
		int [] star;
		ArrayList<Integer> edges = new ArrayList<Integer>();
	}

	public class Edge {
		int v1, v2;
	}

	// use this if you want to store the coordinates
	public boolean store;
	public double[][] coords;
	public double[] fnVertices;
	public int noVertices;
	public int vertexCt;
	public Vertex[] vertices;

	public ArrayList<Triangle> tris = new ArrayList<Triangle>();
	public Triangle [] triangles;
	public int triCt;

	public ArrayList<Edge> edges = new ArrayList<Edge>();
	public ArrayList<Integer> nextTri = new ArrayList<Integer>();
	public ArrayList<Integer> firstTri = new ArrayList<Integer>();
	int edCt;
	int maxStar;
	
	public TriangleDataOpt(boolean storeCoords) {
		store = storeCoords;
	}

	public void addVertex(double xx, double yy, double zz, double ff) {
		if (store) {
			coords[vertexCt][0] = xx;
			coords[vertexCt][1] = yy;
			coords[vertexCt][2] = zz;
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

		// chk if triangle already there
		if (!isAdded(v1, v2, v3)) {
			Triangle t = new Triangle();
			t.v1 = v1;
			t.v2 = v2;
			t.v3 = v3;

			int ein = getPrev(v1, v2);
			if (ein == -1) {
				addNewEdge(v1, v2);
			} else {
				int preTri = nextTri.get(ein);
				nextTri.set(ein, triCt);
				addNextAdj(preTri, v1, v2);
				t.adjTriangle[0][PREV] = preTri;
			}

			ein = getPrev(v2, v3);
			if (ein == -1) {
				addNewEdge(v2, v3);
			} else {
				int preTri = nextTri.get(ein);
				nextTri.set(ein, triCt);
				addNextAdj(preTri, v2, v3);
				t.adjTriangle[1][PREV] = preTri;
			}

			ein = getPrev(v1, v3);
			if (ein == -1) {
				addNewEdge(v1, v3);
			} else {
				int preTri = nextTri.get(ein);
				nextTri.set(ein, triCt);
				addNextAdj(preTri, v1, v3);
				t.adjTriangle[2][PREV] = preTri;
			}
			tris.add(t);
			vertices[v1].starList.add(triCt);
			vertices[v2].starList.add(triCt);
			vertices[v3].starList.add(triCt);
			triCt++;
		}
	}

	private boolean isAdded(int v1, int v2, int v3) {
		int s1 = vertices[v1].starList.size();
		int s2 = vertices[v2].starList.size();
		int s3 = vertices[v3].starList.size();

		if (s1 < s2 && s1 < s3) {
			return isPresent(v1, v1, v2, v3);
		} else if (s2 < s3) {
			return isPresent(v2, v1, v2, v3);
		} else {
			return isPresent(v3, v1, v2, v3);
		}
	}

	private boolean isPresent(int v, int v1, int v2, int v3) {
		int s = vertices[v].starList.size();
		for (int i = 0; i < s; i++) {
			int tin = vertices[v].starList.get(i);
			Triangle t = tris.get(tin);
			if (t.v1 == v1 && t.v2 == v2 && t.v3 == v3) {
				return true;
			}
		}
		return false;
	}

	void addNextAdj(int tin, int v1, int v2) {
		Triangle t = tris.get(tin);
		if (t.v1 == v1) {
			if (t.v2 == v2) {
				// v1 v2
				t.adjTriangle[0][NEXT] = triCt;
			} else {
				// v1 v3
				t.adjTriangle[2][NEXT] = triCt;
			}
		} else {
			// v2 v3
			t.adjTriangle[1][NEXT] = triCt;
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
			Simplex sim = loader.getNextSimplex();
			while (sim != null) {
				if (sim instanceof vgl.iisc.external.types.Vertex) {
					vgl.iisc.external.types.Vertex v = (vgl.iisc.external.types.Vertex) sim;
					double fn;
					int f = Integer.parseInt(ftype);
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
					
				} else if (sim instanceof vgl.iisc.external.types.Edge) {
					pr("Simplicial meshes with edges are not yet supported. This edge will be ignored");
				} else if (sim instanceof vgl.iisc.external.types.Triangle) {
					// TODO Chk if required vertices are added
					vgl.iisc.external.types.Triangle t = (vgl.iisc.external.types.Triangle) sim;
					addTriangle(t.v1, t.v2, t.v3);
				} else if (sim instanceof vgl.iisc.external.types.Tetrahedron) {
					// TODO Chk if required vertices are added
					vgl.iisc.external.types.Tetrahedron t = (vgl.iisc.external.types.Tetrahedron) sim;
					addTetraHedra(t.v1, t.v2, t.v3, t.v4);
				} else {
					er("Invalid Simplex");
				}
				sim = loader.getNextSimplex();
			}
			pr("Finished reading data from file. Loading it......");
			setupTriangles();
			pr("Successfully loaded Data");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	public void setupTriangles() {
		triangles = tris.toArray(new Triangle[0]);
		for(int i = 0;i < edCt;i ++) {
			int tin = nextTri.get(i);
			int fin = firstTri.get(i);
			Edge ein = edges.get(i);
			int ef = 0;
			int et = 0;
			Triangle f = triangles[fin];
			Triangle e = triangles[tin];

			if(ein.v1 == f.v1) {
				if(ein.v2 == f.v2) {
					ef = 0;
				} else {
					ef = 2;
				}
			} else {
				ef = 1;
			}

			if(ein.v1 == e.v1) {
				if(ein.v2 == e.v2) {
					et = 0;
				} else {
					et = 2;
				}
			} else {
				et = 1;
			}

			f.adjTriangle[ef][PREV] = tin;
			e.adjTriangle[et][NEXT] = fin;
		}
		for(int i = 0;i < noVertices;i ++) {
			vertices[i].star = new int[vertices[i].starList.size()];
			int ct = 0;
			for(Iterator<Integer> it = vertices[i].starList.iterator();it.hasNext();) {
				vertices[i].star[ct ++] = it.next(); 
			}
			maxStar = Math.max(maxStar, vertices[i].star.length);
			vertices[i].starList = null;
		}
	}

	public void addTetraHedra(int v1, int v2, int v3, int v4) {
		addTriangle(v1, v2, v3);
		addTriangle(v1, v2, v4);
		addTriangle(v1, v3, v4);
		addTriangle(v2, v3, v4);
	}

	public void setNoOfVertices(int nv) {
		noVertices = nv;
		if (store) {
			coords = new double[noVertices][3];
		}
		fnVertices = new double[noVertices];
		vertices = new Vertex[noVertices];
		
		vertexCt = 0;
		triCt = 0;
		edCt = 0;
	}
	
	public void countTris() {
		int ct = 0;
		for (int i = edCt - 1; i >= 0; i--) {
			int tin = nextTri.get(i);
			Edge ein = edges.get(i);
			do {
				ct++;
				int ed = 0;
				Triangle t = tris.get(tin);
				if (ein.v1 == t.v1) {
					if (ein.v2 == t.v2) {
						ed = 0;
					} else {
						ed = 2;
					}
				} else {
					ed = 1;
				}
				tin = tris.get(tin).adjTriangle[ed][PREV];
			} while (tin != -1);
		}
		System.out.println("No. of tris counted using tris adj. to edges : " + ct / 3);

		ct = 0;
		for (int i = 0; i < vertexCt; i++) {
			ct += vertices[i].starList.size();
		}
		System.out.println("No. of tris counted using star : " + ct / 3);
	}

}
