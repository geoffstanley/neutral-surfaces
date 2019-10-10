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

import java.util.Arrays;

public class MyIntList {
	public int [] array;
	public int length = 0;
	
	public MyIntList(int no) {
		array = new int[no];
	}
	
	public MyIntList() {
		array = new int[10];
	}
	
	public void add(int n) {
		if(array == null) {
			array = new int[10];
		}
		if(length == array.length) {
			array = Arrays.copyOf(array, (int) (length * 1.5));
		}
		array[length ++] = n;
	}
	
	public int size() {
		return length;
	}
	
	public int get(int i) {
		return array[i];
	}

	public void clear() {
		array = null;
		length = 0;
	}

	public void set(int pos, int val) {
		array[pos] = val;		
	}
}
