/*
 *	Copyright (C) 2010 Visualization & Graphics Lab (VGL), Indian Institute of Science
 *
 *	This file is part of libRG, a library to compute Reeb graphs.
 *
 *	libRG is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	libRG is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public License
 *	along with libRG.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Author(s):	Harish Doraiswamy
 *	Version	 :	1.0
 *
 *	Modified by : -- 
 *	Date : --
 *	Changes  : --
 */
package vgl.iisc.external.loader;

import vgl.iisc.external.types.Edge;
import vgl.iisc.external.types.Simplex;
import vgl.iisc.external.types.Tetrahedron;
import vgl.iisc.external.types.Triangle;
import vgl.iisc.external.types.Vertex;

/**
 * 
 * This interface should be extended in order to write a custom loader to support different input file formats.
 *   
 */
public interface MeshLoader {
	
	/**
	 * This method is first called to set the input. You need to perform the required initialization
	 * in this method.
	 * 
	 * @param inputMesh File name of the input mesh
	 */
	public void setInputFile(String inputMesh);
	
	/**
	 * This method is called after the setInputFile() method.
	 * 
	 * @return	The number of vertices in the input. Should be a positive number
	 */
	public int getVertexCount(); 
	
	/**
	 * This method is called until a null is returned, signifying that the mesh is loaded. 
	 * Care should be taken to ensure that the vertices are returned before the edges, 
	 * triangles or tetrahedra incident of that vertex.
	 * 
	 * @return	The next simplex of the input mesh. This can be either {@link Vertex}, {@link Edge}, {@link Triangle} or {@link Tetrahedron}.
	 * 			<br/>
	 * 			null after the last simplex is returned
	 */
	public Simplex getNextSimplex();
	
	/**
	 * This method is called resets the file pointer to the beginning of the input file. 
	 */
	public void reset();

	public int getSimplexCount(); 
	
	public int [] getVertexMap();
}
