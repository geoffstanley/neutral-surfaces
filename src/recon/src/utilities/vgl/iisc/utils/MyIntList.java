/*
 *	Copyright (C) 2012 Visualization & Graphics Lab (VGL), Indian Institute of Science
 *
 *	This file is part of Recon, a library to compute Reeb graphs.
 *
 *	libRG is free software: you can redistribute it and/or modify
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
package vgl.iisc.utils;

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
			array = Arrays.copyOf(array, (int) Math.max(length * 1.5, length + 1));
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
	
	public void remove(int val) {
		int in = -1;
		for(int i = 0;i < length;i ++) {
			if(array[i] == val) {
				in = i;
				break;
			}
		}
		if(in != -1) {
			array[in] = array[length - 1];
			length --;
		}
	}
	
	public int removeIndex(int index) {
		int in = index;
		int ret = -1;
		if(in != -1 && in >= 0 && in < length) {
			ret = array[in];
			array[in] = array[length - 1];
			length --;
		}
		return ret;
	}
	
	public void removeElement(int el) {
		int i = 0;
		for(;i < length;i ++) {
			if(array[i] == el) {
				break;
			}
		}
		if(i == length) {
			return;
		}
		for(;i < length - 1;i ++) {
			array[i] = array[i + 1];
		}
		length --;
	}
	
	public void replace(int replacee, int replacer) {
		for(int i = 0;i < length;i ++) {
			if(array[i] == replacee) {
				array[i] = replacer;
				return;
			}
		}
		System.out.println("Could not find element!!");		
	}
	
	public boolean contains(int el) {
		for(int i = 0;i < length;i ++) {
			if(array[i] == el) {
				return true;
			}
		}
		return false;
	}
	
	public void addAll(MyIntList list) {
		for(int i = 0;i < list.length;i ++) {
			add(list.get(i));
		}
	}
}
