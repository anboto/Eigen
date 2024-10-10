// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _Eigen_Eigen_h
#define _Eigen_Eigen_h

#define EIGEN_MATRIX_PLUGIN 	<Eigen/ToStringPlugin.h>
#define EIGEN_DENSEBASE_PLUGIN 	<Eigen/ToStringPlugin.h>
#define EIGEN_TENSOR_PLUGIN		<Eigen/ToStringPlugin.h>

#define EIGEN_MPL2_ONLY

#ifndef _DEBUG
#define EIGEN_NO_DEBUG
#else
#define EIGEN_INITIALIZE_MATRICES_BY_NAN 
#endif

#define eigen_assert(x) ASSERT(x)

#undef Success  
#include <plugin/eigen/Eigen/Dense>
#include <plugin/eigen/unsupported/Eigen/NonLinearOptimization>
#undef Complex
#include <plugin/eigen/unsupported/Eigen/FFT>
#include <plugin/eigen/unsupported/Eigen/CXX11/Tensor>

#include "MultiDimMatrix.h"

namespace Upp {

template <class T>
using UVector = Upp::Vector<T>;

template <class T>
using UArray = Upp::Array<T>;

template <class T>
using UIndex = Upp::Index<T>;

template<typename _Scalar, ptrdiff_t nx = Eigen::Dynamic, ptrdiff_t ny = Eigen::Dynamic>
struct NonLinearOptimizationFunctor {
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = nx,
		ValuesAtCompileTime = ny
	};
	typedef Eigen::Matrix<double, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix<double, ValuesAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix<double, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;
	
	Eigen::Index unknowns, datasetLen;
	
	NonLinearOptimizationFunctor() : unknowns(InputsAtCompileTime), datasetLen(ValuesAtCompileTime) {}
	NonLinearOptimizationFunctor(int unknowns, int datasetLen) : unknowns(unknowns), datasetLen(datasetLen) {}
	
	ptrdiff_t inputs() const {return ptrdiff_t(unknowns);}
	ptrdiff_t values() const {return ptrdiff_t(datasetLen);}
	virtual void operator() (const InputType& , ValueType* , JacobianType*  = 0) const {};
};

struct Basic_functor : NonLinearOptimizationFunctor<double> {
	Basic_functor(Function <int(const Eigen::VectorXd &b, Eigen::VectorXd &err)> _function) : function(_function) {}
	int operator()(const Eigen::VectorXd &b, Eigen::VectorXd &fvec) const {return function(b, fvec);}
	Function <int(const Eigen::VectorXd &b, Eigen::VectorXd &err)> function;
};

bool NonLinearOptimization(Eigen::VectorXd &y, Eigen::Index numData, 
			Function <int(const Eigen::VectorXd &y, Eigen::VectorXd &residual)>residual,
			double xtol = Null, double ftol = Null, int maxfev = Null);
bool SolveNonLinearEquations(Eigen::VectorXd &y, Function <int(const Eigen::VectorXd &b, Eigen::VectorXd &residual)> Residual,
			double xtol = Null, int maxfev = Null, double factor = Null);
double SolveNonLinearEquation(double y, Function <double(double b)> Residual, double xtol = Null, int maxfev = Null, double factor = Null);

Eigen::Matrix3d SkewSymmetricMatrix(const Eigen::Vector3d& r);		// Antisymmetric matrix from 3d vector
	
template <class T>
void Xmlize(XmlIO &xml, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Size_<int64> sz(mat.cols(), mat.rows());
	xml ("size", sz);
	if(xml.IsStoring()) {
		for(int r = 0; r < mat.rows(); r++)
			for(int c = 0; c < mat.cols(); c++) {
				XmlIO io = xml.Add("item");
				T data = mat(r, c);
				Xmlize(io, data);
			}
	} else {
		mat.resize(ptrdiff_t(sz.cy), ptrdiff_t(sz.cx));
		int r = 0, c = 0;
		for(int i = 0; i < xml->GetCount(); i++) 
			if(xml->Node(i).IsTag("item")) {
				XmlIO io = xml.At(i);
				T data;
				Xmlize(io, data);
				mat(r, c) = data;
				++c;
				if (c == sz.cx) {
					c = 0;
					r++;
				}
			}
	}
}

template <class T>
void Xmlize(XmlIO &xml, Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
	int64 sz = vec.size();
	xml ("size", sz);
	if(xml.IsStoring()) {
		for(int r = 0; r < sz; r++) {
			XmlIO io = xml.Add("item");
			T data = vec(r);
			Xmlize(io, data);
		}
	} else {
		vec.resize(ptrdiff_t(sz));
		int r = 0;
		for(int i = 0; i < xml->GetCount(); i++)
			if(xml->Node(i).IsTag("item")) {
				XmlIO io = xml.At(i);
				T data;
				Xmlize(io, data);
				vec(r++) = data;
			}
	}
}

template <typename T, int NumIndices>
void Jsonize(JsonIO &io, Eigen::Tensor<T, NumIndices> &mat) {
	Array<T> vector;
	Vector<int64> vsz(NumIndices);
	for (int i = 0; i < NumIndices; ++i)
		vsz[i] = mat.dimension(i);
	io("size", vsz);
	if(io.IsStoring()) {
		vector.SetCount(int(mat.size()));
		Copy(mat, vector);
		io("data", vector);
	} else {
		io("data", vector);
		Eigen::array<Eigen::Index, NumIndices> dims;
		for (int i = 0; i < NumIndices; ++i) 	
			dims[i] = static_cast<Eigen::Index>(vsz[i]);
		mat.resize(dims);
		std::copy(vector.begin(), vector.end(), mat.data());
	}
}

template <class T>
void Jsonize(JsonIO &io, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Array<T> vector;
	Size_<int64> sz(mat.cols(), mat.rows());
	io("size", sz);
	if(io.IsStoring()) {
		vector.SetCount(int(sz.cx)*int(sz.cy));
		Copy(mat, vector);
		io("data", vector);
	} else {
		io("data", vector);
		mat.resize(ptrdiff_t(sz.cy), ptrdiff_t(sz.cx));
		std::copy(vector.begin(), vector.end(), mat.data());
	}
}

template <class T>
void Jsonize(JsonIO &io, Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
	Array<T> vector;
	int64 sz = vec.size();
	io("size", sz);
	if(io.IsStoring()) {
		vector.SetCount(int(sz));
		Copy(vec, vector);
		io("data", vector);
	} else {
		io("data", vector);
		vec.resize(ptrdiff_t(sz));
		Copy(vector, vec);
	}
}

template <class T>
void Serialize(Stream& stream, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Size_<int64> sz(mat.cols(), mat.rows());
	stream % sz;
	if(stream.IsStoring()) {
		for(int r = 0; r < mat.rows(); r++)
			for(int c = 0; c < mat.cols(); c++) {
				T data = mat(r, c);
				stream % data;
			}
	} else {
		mat.resize(ptrdiff_t(sz.cy), ptrdiff_t(sz.cx));
		int r = 0, c = 0;
		for(int i = 0; i < sz.cy*sz.cx; i++) {
			T data;
			stream % data;
			mat(r, c) = data;
			++c;
			if (c == sz.cx) {
				c = 0;
				r++;
			}
			if (r == sz.cy)
				break;
		}
	}
}

template <class T>
void Serialize(Stream& stream, Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
	int64 sz = vec.size();
	stream % sz;
	if(stream.IsStoring()) {
		for (int i = 0; i < sz; ++i) {
			T data = vec(i);
			stream % data;
		}
	} else {
		vec.resize(ptrdiff_t(sz));
		for (int i = 0; i < sz; ++i) {
			T data;
			stream % data;
			vec(i) = data;
		}
	}
}

// These functions serve both for Eigen, std and U++ Vectors

template <class Range>
void Resize(Range &v, size_t len) {v.SetCount(int(len));}
template <class Range>
void Resize(Range &v, size_t len, const typename Range::value_type& init) {
	v.SetCount(int(len));
	std::fill(v.begin(), v.end(), init);
}

template <class Range>
void ResizeConservative(Range &v, size_t len) {v.SetCount(int(len));}
template <class Range>
void ResizeConservative(Range &v, size_t len, const typename Range::value_type& init) {v.SetCount(int(len), init);}
template <class Range>
void Clear(Range &v) {v.Clear();}

template <typename T>
void Resize(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len) {v.resize(len);}
template <typename T>
void Resize(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len, const T& init) {v.setConstant(len, 1, init);}
template <typename T>
void ResizeConservative(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len) {v.conservativeResize(len);}
template <typename T>
void ResizeConservative(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len, const T& init) {
	size_t len0 = v.size();
	v.conservativeResize(len);
	if (len > len0)
		std::fill(&v[len0], v.data() + len, init);
}
template <typename T>
void Clear(Eigen::Matrix<T, Eigen::Dynamic, 1> &v) 				{v.resize(0);}
template <typename T>
void Clear(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &m) {m.resize(0, 0);}
template <typename T, int NumIndices_>
void Clear(Eigen::Tensor<T, NumIndices_> &m) 					{m = Eigen::Tensor<T, NumIndices_>();}

template <typename T>
void PrePad(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, size_t len, const T& init) {
	size_t len0 = v.size();
	v.conservativeResize(len);
	if (len > len0) {
		size_t delta = len - len0;
		std::copy(&v[len0 - delta], v.data() + len0, &v[len0]);
		std::copy(v.data(), v.data() + len0 - delta, &v[delta]);
		std::fill(v.data(), v.data() + delta, init);
	}
}


template <typename T>
void Resize(std::vector<T> &v, size_t len) {v.resize(len);}
template <typename T>
void Resize(std::vector<T> &v, size_t len, const T& init) {
	v.resize(len);
	std::fill(v.begin(), v.end(), init);
}
template <typename T>
void ResizeConservative(std::vector<T> &v, size_t len) {v.resize(len);}
template <typename T>
void ResizeConservative(std::vector<T> &v, size_t len, const T& init) {v.resize(len, init);}
template <typename T>
void Clear(std::vector<T> &v) {v.clear();}

#define PostPad ResizeConservative


template <typename T>
auto Begin(const std::vector<T> &v)		{return v.begin();}
template <typename T>
auto Begin(std::vector<T> &v)			{return v.begin();}
template <typename T>
auto End(const std::vector<T> &v)		{return v.end();}
template <typename T>
auto End(std::vector<T> &v)				{return v.end();}

template <typename T>
auto Begin(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v){return v.data();}
template <typename T>
auto Begin(Eigen::Matrix<T, Eigen::Dynamic, 1> &v)		{return v.data();}
template <typename T>
auto End(const Eigen::Matrix<T, Eigen::Dynamic, 1> &v)	{return v.data() + v.size();}

template <typename T>
auto Begin(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v)	{return v.data();}
template <typename T>
auto Begin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v)			{return v.data();}
template <typename T>
auto End(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v)		{return v.data() + v.size();}

template <typename T, int NumIndices>
auto Begin(const Eigen::Tensor<T, NumIndices> &v)		{return v.data();}
template <typename T, int NumIndices>
auto Begin(Eigen::Tensor<T, NumIndices> &v)				{return v.data();}
template <typename T, int NumIndices>
auto End(const Eigen::Tensor<T, NumIndices> &v)			{return v.data() + v.size();}

template <class Range>
auto Begin(const Range &v)				{return v.Begin();}
template <class Range>
auto Begin(Range &v)					{return v.Begin();}
template <class Range>
auto End(const Range &v)				{return v.End();}
template <class Range>
auto End(Range &v)						{return v.End();}


template <class Range>
auto &First(Range &data) {return data[0];}

template <class Range>
auto &Last(Range &data) {return data[data.size()-1];}

template <typename T>			// To avoid problem with Upp
void ReverseX(std::vector<T> &v) {std::reverse(v.begin(), v.end());}

template <typename T>			// To avoid problem with Upp
void ReverseX(Eigen::Matrix<T, Eigen::Dynamic, 1> &v) {v.reverseInPlace();}

template <typename T>			// To avoid problem with Upp
void ReverseX(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &v) {v.reverseInPlace();}

template <class Range>
void ReverseX(Range &v) {		// To avoid problem with Upp
	typename Range::value_type *first = v.begin();
	typename Range::value_type *last = v.end();
	while ((first != last) && (first != --last)) 
		Swap(*first++, *last);
}

template <typename T>			
void Rotate(std::vector<T> &v, int shift) {
	if (shift > 0)
		std::rotate(v.begin(), v.begin() + shift, v.end());
	else if (shift < 0)
		std::rotate(v.rbegin(), v.rbegin() - shift, v.rend());
}

template <class Range>			
void Rotate(Range &v, int k) {
	auto rotate = [](Range &v, int start, int end) {
	    while (start < end) {
	        Swap(v[start], v[end]);
	        start++;
	        end--;
	    }
	};
	int n = v.size();
	k = k % n;  // Handle cases where k is greater than the array size
	if (k > 0) {
	    rotate(v, 0, n - 1);
	    rotate(v, 0, k - 1);
	    rotate(v, k, n - 1);
	} else if (k > 0) {
	    rotate(v, 0, k - 1);
	    rotate(v, k, n - 1);
	    rotate(v, 0, n - 1);
	}
}

template <class Range>
void CopyRowMajor(Range &in, int nrows, int ncols, Eigen::Matrix<typename Range::value_type, Eigen::Dynamic, Eigen::Dynamic> &out) {
	out = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(Begin(in), nrows, ncols);
}

template <typename T>
void CopyRowMajor(T *in, int nrows, int ncols, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &out) {
	out = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(in, nrows, ncols);
}

template <class Range>
void CopyRowMajor(const Eigen::Matrix<typename Range::value_type, Eigen::Dynamic, Eigen::Dynamic> &in, Range &out) {
	Resize(out, in.size());
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajMat;
	RowMajMat::Map(Begin(out), in.rows(), in.cols()) = in;
}

template <class Range1, class Range2>
void Copy(const Range1& in, Range2 &out) {
	Resize(out, in.size());
	std::copy(Begin(in), End(in), Begin(out));
}

template <class Range1, class Range2>
void Block(Range1& in, Range2& out, int begin, int len) {
	ASSERT(begin < in.size() && len > 0 && begin + len <= in.size());
	Resize(out, len);
	std::copy(in + begin, in + begin + len, Begin(out));
}


template <class T>
void Swap(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, int rc1, int rc2) {
 	A.row(rc1).swap(A.row(rc2));	// Swap rows rc1 and rc2
 	A.col(rc1).swap(A.col(rc2));	// Swap columns rc1 and rc2
}

template <class T>
void Swap(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A1, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A2, int rc1, int rc2) {
	ASSERT(A1.rows() == A2.rows() && A1.cols() == A2.cols());
	
	A1.row(rc1).swap(A2.row(rc2));
    A1.col(rc1).swap(A2.col(rc2));	
}

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> Segment(const Eigen::Matrix<T, Eigen::Dynamic, 1> &d, int ifrom, int num) {
	return d.segment(ifrom, num);
}

template <typename T>
inline std::vector<T> Segment(const std::vector<T> &d, int ifrom, int num) {
	return std::vector<T>(d.begin() + ifrom, d.begin() + ifrom + num);
}

template <class Range>
inline Range Segment(const Range &d, int ifrom, int num) {
	Range a;
	if (ifrom + num >= d.size()) {
		num = d.size() - ifrom;
		if (num <= 0)
			return a;
	}
	Resize(a, num);
	std::copy(Begin(d) + ifrom, Begin(d) + ifrom + num, Begin(a));
	return a;
}

template <class Range>
inline Range Mid(const Range &d, int ifrom, int num) {return Segment(d, ifrom, num);}

template <class Range>
inline Range Left(const Range &d, int num) {
	Range a;
	if (num >= d.size()) {
		a = clone(d);
		return a; 
	}
	Resize(a, num);
	std::copy(Begin(d), Begin(d) + num, Begin(a));
	return a;
}

template <class Range>
inline Range Right(const Range &d, int num) {
	Range a;
	if (num >= d.size()) {
		a = clone(d);
		return a; 
	}
	Resize(a, num);
	std::copy(Begin(d) + d.size() - num, Begin(d) + d.size(), Begin(a));
	return a;
}

template <typename T>
inline void Remove(Vector<T> &d, int id) {d.Remove(id);}

template <typename T>
inline void Remove(Eigen::Matrix<T, Eigen::Dynamic, 1> &d, int id) {
	Eigen::Index sz = d.size();
	Eigen::Index right = sz-id-1;
	d.segment(id, right) = d.segment(id+1, right);
	d.conservativeResize(sz-1);
}

template <class Range>
inline void Remove(Range &d, int id) {
	d.erase(d.begin() + id);
}

template <class T>
bool IsNull(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &a) {return a.size() == 0;}

#define EigenNull	Eigen::MatrixXd()

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;
   
template<typename Scalar, int rank, typename sizeType>
auto TensorToMatrix(const Eigen::Tensor<Scalar,rank> &tensor, const sizeType rows, const sizeType cols) { 
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows, cols);
}

template<typename Scalar, typename... Dims>
auto MatrixToTensor(const MatrixType<Scalar> &matrix, Dims... dims) {
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}

template<typename T>
decltype(auto) TensorLayoutSwap(T&& t) {
	return Eigen::TensorLayoutSwapOp<typename std::remove_reference<T>::type>(t);
}
 
}


namespace Eigen {

template<> struct NumTraits<Upp::Complex> : GenericNumTraits<Upp::Complex>
{
  typedef double Real;
  typedef typename NumTraits<double>::Literal Literal;
  enum {
    IsComplex = 1,
    RequireInitialization = NumTraits<Real>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real epsilon() { return NumTraits<Real>::epsilon(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR
  static inline int digits10() { return NumTraits<Real>::digits10(); }
};

}
#endif