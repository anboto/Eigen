// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _ScatterDraw_MultiDimMatrixIndex_h_
#define _ScatterDraw_MultiDimMatrixIndex_h_

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
	int GetNumData() const {
		int ret = 1;
		for (auto dim : axisDim)
			ret *= dim;
		return ret;
	}
	int size() const			{return GetNumData();}
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
	
	inline int GetNumAxis() const	{return index.GetNumAxis();}
		
    T& operator()(int row, int col) {
        ASSERT(index.IsValid(row, col));
        return d[index(row, col)];
    }
    const T& operator()(int row, int col) const {
        ASSERT(index.IsValid(row, col));
        return d[index(row, col)];
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
	
}

#endif
