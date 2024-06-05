// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2024, the Anboto author and contributors
#ifndef _ScatterDraw_MultiDimMatrix_h_
#define _ScatterDraw_MultiDimMatrix_h_

#include <Eigen/Eigen.h>

namespace Upp {

class MultiDimMatrixIndex {
public:
	MultiDimMatrixIndex()			   {};
	template<typename... Args>
	MultiDimMatrixIndex(Args... args)  {SetAxis(args...);}
	
	void SetNumAxis(int numAxis)	   {axisDim.SetCount(numAxis);};
	inline int GetNumAxis() const  	   {return axisDim.size();}
	
	void SetAxisDim(int axis, int dim) {
		ASSERT(axis >= 0 && axis < axisDim.size() && dim > 0);
		axisDim[axis] = dim;
	}
	void InsertAxis(int axis, int dim) {
		axisDim.Insert(axis, dim);
	}
	template<typename... Args>
	MultiDimMatrixIndex &SetAxis(int t, Args... args) {
		ASSERT(t > 0);
		axisDim << t;
		SetAxis(args...);
		return *this;
	}
	MultiDimMatrixIndex &SetAxis(const Vector<int> &dim) {
		axisDim = clone(dim);
		return *this;
	}
			
	Vector<int> &GetAxisDim()	{return axisDim;};
	
	int GetIndex(const Vector<int> &idx) const {
		ASSERT(IsValid(idx));
		if (colMajor) {
			int index = 0;
			int multiplier = 1;
			for (int ix = 0; ix < axisDim.size(); ++ix) {
				index += multiplier*idx[ix];
				multiplier *= axisDim[ix];
			}
			return index;
		} else {
			int index = 0;
			int multiplier = 1;
			for (int ix = axisDim.size()-1; ix >= 0; --ix) {
				index += multiplier*idx[ix];
				multiplier *= axisDim[ix];
			}
			return index;
		}
	}
	template<typename T, typename... Args>
	int GetIndex(T t, Args... args) const {
		Vector<int> index;
		
		index << t;
	    AddIndex(index, args...);
	
	    return GetIndex(index);
	}
	inline int GetIndex(int row, int col) const {
		ASSERT(IsValid(row, col));
		if (colMajor) 
			return row + axisDim[0]*col;
		else
			return col + axisDim[1]*row;
	}
	
	void ResetIndex(Vector<int> &idx) const {
		idx.SetCount(axisDim.size(), 0);
	}
	bool IncrementIndex(Vector<int> &idx) const {
		ASSERT(idx.size() == axisDim.size());
		for (int i = 0; i < axisDim.size(); ++i) {
			idx[i]++;
			if (idx[i] < axisDim[i])
				return true;
			idx[i] = 0;
		}
		return false;
	}
		
	template<typename... Args>
	inline int operator()(Args... args) const 	{return GetIndex(args...);}
	
	inline int operator()(int row, int col) const  	{return GetIndex(row, col);}
		
	bool IsValid(const Vector<int> &idx) const {		
		for (int ix = 0; ix < axisDim.size(); ++ix) 
			if (idx[ix] < 0 && idx[ix] >= axisDim[ix])
				return false;
		return true;
	}
	template<typename T, typename... Args>
	bool IsValid(T t, Args... args) const {
		Vector<int> index;
		
		index << t;
	    AddIndex(index, args...);
	
	    return IsValid(index);
	}	
	inline bool IsValid(int row, int col) const  {
		return row >= 0 && row < axisDim[0] && col >= 0 && col < axisDim[1];
	}
	int size() const {
		int ret = 0;
		for (auto dim : axisDim)
			if (dim > 0) {
				if (ret == 0)
					ret = 1;
				ret *= dim;
			}
		return ret;
	}
	int size(int dim) const		{return axisDim[dim];}
	
	MultiDimMatrixIndex &ColMajor(bool c = true)	{colMajor = c;	return *this;}
	MultiDimMatrixIndex &RowMajor(bool c = true)	{colMajor = !c;	return *this;}
	
	void Xmlize(XmlIO xml) {
		xml
			("axisDim", axisDim)	
		;	
	}
	void Jsonize(JsonIO& json) {
		json
			("axisDim", axisDim)	
		;	
	}
private:
	Vector<int> axisDim;
	
	template<typename T>
	static void AddIndex(Vector<int> &index, T t) {
		index << t;
	}	
	template<typename T, typename... Args>
	static void AddIndex(Vector<int> &index, T t, Args... args) {
		index << t;
		AddIndex(index, args...);
	}
	void SetAxis(int dimX) {
		ASSERT(dimX > 0);
		axisDim << dimX;
	}

protected:
	bool colMajor = true;
};

class MultiDimMatrixIndexRowMajor : public MultiDimMatrixIndex {
public:
	MultiDimMatrixIndexRowMajor() 				{colMajor = false;};
	template<typename... Args>
	MultiDimMatrixIndexRowMajor(Args... args) 	{colMajor = false; SetAxis(args...);}
};

template <class T>
class MultiDimMatrix {
public:
	MultiDimMatrix()			  	{};
	template<typename... Args>
	MultiDimMatrix(Args... args)  	{index.SetAxis(args...);}
	
	template<typename... Args>
	void Resize(int t, Args... args) {
		index.SetAxis(t, args...);
		d.Alloc(index.size());
	}
	
	void Resize(const Vector<int> &dim) {
		index.SetAxis(dim);
		d.Alloc(index.size());
	}
	
	void InsertAxis(int axis, int dim) {
		index.InsertAxis(axis, dim);
	}
	
	inline int GetNumAxis() const	{return index.GetNumAxis();}
		
    T& operator()(int row, int col) {
        ASSERT(index.IsValid(row, col));
        return d[index(row, col)];
    }
    const T& operator()(int row, int col) const {
        ASSERT(index.IsValid(row, col));
        return d[index(row, col)];
    }
    
    T& operator()(const Vector<int> &indx) {
        ASSERT(index.IsValid(indx));
        return d[index.GetIndex(indx)];
    }
    const T& operator()(const Vector<int> &indx) const {
        ASSERT(index.IsValid(indx));
        return d[index.GetIndex(indx)];
    }
    
    template<typename... Args>
	T &operator()(Args... args) {
		ASSERT(index.IsValid(args...));
		return d[index.GetIndex(args...)];
	}
	template<typename... Args>
	const T &operator()(Args... args) const {
		ASSERT(index.IsValid(args...));
		return d[index.GetIndex(args...)];
	}
	
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> block(const Vector<int> &indx, int idrow, int wrow, int idcol, int wcol) {
		ASSERT(indx.size() == GetNumAxis());
		ASSERT(idrow >= 0 && idcol >= 0);
		ASSERT(indx[idrow]+wrow < size(idrow));
		ASSERT(indx[idcol]+wcol < size(idcol));
		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ret(wrow, wcol);	
		Vector<int> id = clone(indx);
		
		for (int r = 0; r < wrow; ++r, id[idrow]++) {
			id[idcol] = indx[idcol];
			for (int c = 0; c < wcol; ++c, id[idcol]++) 
				ret(r, c) = d[index.GetIndex(id)];
		}
		return ret;
	}
	
	const T *begin() const		{return d.begin();}
	T *begin() 					{return d.begin();}
	
	void ColMajor(bool c = true){index.ColMajor(c);}
	void RowMajor(bool c = true){index.RowMajor(c);}
	
	int size() const			{return index.size();}
	int size(int dim) const		{return index.size(dim);}

protected:
    Buffer<T> d;
	MultiDimMatrixIndex index;
};

template <class T>
class MultiDimMatrixRowMajor : public MultiDimMatrix<T> {
public:
	MultiDimMatrixRowMajor() 				{this->RowMajor();};
	template<typename... Args>
	MultiDimMatrixRowMajor(Args... args) 	{this->RowMajor(); this->Resize(args...);}
};

template <typename T>
void RowMajorToColMajor(const T *d_row, T *d_col, const Vector<int> &dims) {
	MultiDimMatrixIndex row, col;
	row.RowMajor().SetAxis(dims);
	col.ColMajor().SetAxis(dims);
	Vector<int> idx;
	row.ResetIndex(idx);
	do {
		d_col[col.GetIndex(idx)] = d_row[row.GetIndex(idx)];
	} while(row.IncrementIndex(idx));
}

template <typename T>
void ColMajorToRowMajor(const T *d_col, T *d_row, const Vector<int> &dims) {
	MultiDimMatrixIndex row, col;
	row.RowMajor().SetAxis(dims);
	col.ColMajor().SetAxis(dims);
	Vector<int> idx;
	row.ResetIndex(idx);
	do {
		d_row[row.GetIndex(idx)] = d_col[col.GetIndex(idx)];
	} while(row.IncrementIndex(idx));
}

}

#endif
