// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2025, the Anboto author and contributors
#ifndef _ScatterDraw_MultiDimMatrix_h_
#define _ScatterDraw_MultiDimMatrix_h_

#include <Eigen/Eigen.h>

namespace Upp {

class MultiDimMatrixIndex {
public:
	MultiDimMatrixIndex()			   	{};
	template<typename... Args>
	MultiDimMatrixIndex(Args... args)  	{SetAxis(args...);}
	
	void SetNumAxis(int numAxis)	   	{axisDim.SetCount(numAxis);};
	inline int GetNumAxis() const  	   	{return axisDim.size();}
	const Vector<int> &GetAxisDim()const{return axisDim;}
	int GetAxisDim(int dim)	const		{return axisDim[dim];}
	
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
	template<typename T>
	int GetIndex(T t) const {
		Vector<int> index;
		
	    AddIndex(index, t);
	
	    return GetIndex(index);
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
	Vector<int> ResetIndex() const {
		Vector<int> idx;
		idx.SetCount(axisDim.size(), 0);
		return idx;
	}
	void Clear()	{axisDim.Clear();}
	
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
	
	bool IncrementIndexRows(Vector<int> &idx) const {
		return IncrementIndexRowsDim(idx) >= 0;
	}	
	int IncrementIndexRowsDim(Vector<int> &idx) const {
		ASSERT(idx.size() == axisDim.size());
		for (int i = axisDim.size()-1; i >= 0; --i) {
			idx[i]++;
			if (idx[i] < axisDim[i])
				return i;
			idx[i] = 0;
		}
		return -1;
	}
		
	template<typename... Args>
	inline int operator()(Args... args) const 		{return GetIndex(args...);}
	
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
	inline bool IsValid(int row) const  {
		return row >= 0 && row < axisDim[0];
	}
	int size() const {
		if (IsEmpty())
			return 0;
		int ret = 1;
		for (auto dim : axisDim)
			ret *= dim;
		return ret;
	}
	int size(int dim) const		{return axisDim[dim];}
	bool IsEmpty() const		{return axisDim.IsEmpty();}
	
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
	MultiDimMatrix(Args... args) {
		index.SetAxis(args...);
		d.Alloc(index.size());
	}
	
	template<typename... Args>
	void Resize(int t, Args... args) {
		index.SetAxis(t, args...);
		d.Alloc(index.size());
	}
	
	void Resize(const Vector<int> &dim) {
		index.SetAxis(dim);
		d.Alloc(index.size());
	}
	
	void Clear() {
		index.Clear();
		d.Clear();
	}
	
	void InsertAxis(int axis, int dim) {
		index.InsertAxis(axis, dim);
	}
	
	void SetNumAxis(int numAxis)	   	{index.SetNumAxis(numAxis);};
	inline int GetNumAxis() const  	   	{return index.GetNumAxis();}
	const Vector<int> &GetAxisDim()const{return index.GetAxisDim();}
	int GetAxisDim(int dim)	const		{return index.GetAxisDim(dim);}
	
	MultiDimMatrix &SetZero() {
		memset(d.begin(), 0, size()*sizeof(T));
		return *this;
	}
	MultiDimMatrix &SetConstant(const T &val) {
		int sz = size();
		for (int i = 0; i < sz; ++i)
			d[i] = val;
		return *this;
	}
	
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
	
	String ToString() const {
		String r;
	    Vector<int> idx = index.ResetIndex();
	    int num_dims = index.GetNumAxis();
	
	    Function<void(int, int)> PrintRecursively;
	    PrintRecursively = [&](int dim, int elem_count) {
	        if (dim > 0 && r[r.GetCount()-1] != '[') 
	        	r << "\n" << String(' ', dim);
	        r << "[";
	
	        if (dim == num_dims - 1) {
	            for (int i = 0; i < index.size(dim); i++) {
	                if (i) 
	                	r << " ";
	                r << d[elem_count + i];
	            }
	        } else {
	            for (int i = 0; i < index.size(dim); i++) 
	                PrintRecursively(dim + 1, elem_count + i * index.size(dim + 1));
	        }
	        r << "]";
	    };
	    PrintRecursively(0, 0);
	    return r;
	}
		
	const T *begin() const		{return d.begin();}
	T *begin() 					{return d.begin();}
	
	const T &array(int i) const	{return d[i];}
	T &array(int i)				{return d[i];}
	
	void ColMajor(bool c = true){index.ColMajor(c);}
	void RowMajor(bool c = true){index.RowMajor(c);}
	
	int size() const			{return index.size();}
	int size(int dim) const		{return index.size(dim);}
	bool IsEmpty() const		{return index.IsEmpty();}
	
	bool IsNullInstance() const;

	void Xmlize(XmlIO xml) {
		xml
			("index", index)
			("d", d)	
		;	
	}
	void Jsonize(JsonIO& json) {
		json
			("index", index)
			("d", d)	
		;	
	}
protected:
    Buffer<T> d;
	MultiDimMatrixIndex index;
};

template <>
bool MultiDimMatrix<Eigen::Matrix<double, -1, -1>>::IsNullInstance() const;

template <>
bool MultiDimMatrix<Eigen::Matrix<double, -1, 1>>::IsNullInstance() const;


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
	Vector<int> idx = row.ResetIndex();
	do {
		d_col[col.GetIndex(idx)] = d_row[row.GetIndex(idx)];
	} while(row.IncrementIndex(idx));
}

template <typename T>
void ColMajorToRowMajor(const T *d_col, T *d_row, const Vector<int> &dims) {
	MultiDimMatrixIndex row, col;
	row.RowMajor().SetAxis(dims);
	col.ColMajor().SetAxis(dims);
	Vector<int> idx = row.ResetIndex();
	do {
		d_row[row.GetIndex(idx)] = d_col[col.GetIndex(idx)];
	} while(row.IncrementIndex(idx));
}

template <typename T>
bool IsNum(const MultiDimMatrix<T> &a) {
	return !a.IsNullInstance();
}

}

#endif
