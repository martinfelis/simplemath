#pragma once

#include <sstream>
#include <cassert>
#include <iostream>
#include <cstring> // for memcpy

namespace SimpleMath {

template <typename Derived, typename ScalarType, int NumRows, int NumCols>
struct Block;

template <typename ScalarType, int NumRows, int NumCols>
struct Dynamic;

template <typename Derived, typename ScalarType, int NumRows, int NumCols>
struct Transpose;

template <typename Derived, typename ScalarType, int Rows, int Cols>
struct MatrixBase {
	template <typename OtherDerived>
	Derived& operator=(const OtherDerived& other) {
		Derived result (other.rows(), other.cols());

		unsigned int i,j;
		for (i = 0; i < rows(); i++) {
			for (j = 0; j < cols(); j++) {
				this->operator()(i,j) = other(i,j);
			}
		}
		
		return result;
	}
	
	template <typename OtherDerived>
	bool operator==(const OtherDerived& other) {
		unsigned int i,j;
		for (i = 0; i < rows(); i++) {
			for (j = 0; j < cols(); j++) {
				if (this->operator()(i,j) != other(i,j))
					return false;
			}
		}
		
		return true;
	}
	
	template <typename OtherDerived>
	Derived operator+(const OtherDerived& other) {
		Derived result (*(static_cast<Derived*>(this)));
		result += other;
		return result;
	}

	template <typename OtherDerived>
	Derived operator*(const OtherDerived& other) {
		Derived result (Derived::Zero(rows(), other.cols()));

		unsigned int i,j,k;
		for (i = 0; i < rows(); i++) {
			for (j = 0; j < other.cols(); j++) {
				for (k = 0; k < other.rows(); k++) {
					result (i,j) += operator()(i,k) * other(k,j);
				}
			}
		}
		return result;
	}

	size_t rows() const {
		return static_cast<const Derived*>(this)->rows();
	}
	size_t cols() const {
		return static_cast<const Derived*>(this)->cols();
	}

	const ScalarType& operator()(const size_t& i, const size_t& j) const {
		return static_cast<const Derived*>(this)->operator()(i,j);
	}
	ScalarType& operator()(const size_t& i, const size_t& j) {
		return static_cast<Derived*>(this)->operator()(i,j);
	}

	ScalarType* data() {
		return static_cast<Derived*>(this)->data();
	}

	const ScalarType* data() const {
		return static_cast<const Derived*>(this)->data();
	}

	template <
		int block_rows, 
		int block_cols
		>
	Block<
		Derived,
		ScalarType,
		block_rows,
		block_cols
		> block(int block_row_index, int block_col_index) {
			assert(block_row_index + block_rows <= rows());
			assert(block_col_index + block_cols <= cols());
		return Block<Derived, ScalarType, block_rows, block_cols>(static_cast<Derived*>(this), block_row_index, block_col_index);
	}

	Transpose<Derived, ScalarType, Rows, Cols> transpose() {
		return Transpose<Derived, ScalarType, Rows, Cols>(static_cast<Derived*>(this));
	}

	template <typename OtherDerived>
	Transpose<OtherDerived, ScalarType, Rows, Cols> transpose() {
		return Transpose<OtherDerived, ScalarType, Rows, Cols>(static_cast<Derived*>(this));
	}


	static Derived Zero(int NumRows, int NumCols) {
		Derived result (NumRows, NumCols);

		for (size_t i = 0; i < NumRows; i++) {
			for (size_t j = 0; j < NumCols; j++) {
				result(i,j) = static_cast<ScalarType>(0.0);
			}
		}

		return result;
	}
};

template <typename Derived, typename ScalarType, int NumRows, int NumCols>
struct Transpose : public MatrixBase<Transpose<Derived, ScalarType, NumRows, NumCols>, ScalarType, NumRows, NumCols> {
	typedef MatrixBase<Derived, ScalarType, NumRows, NumCols> MatrixType;
	Derived* mTransposeSource;	

	Transpose(Derived* transpose_source) :
		mTransposeSource(transpose_source)
	{ }

	Transpose(const Transpose &other) :
		mTransposeSource(other.mTransposeSource)
	{ }

	template <typename OtherDerived>
	Transpose& operator=(const OtherDerived& other) {
		unsigned int i,j;
		for (i = 0; i < rows(); i++) {
			for (j = 0; j < cols(); j++) {
				this->operator()(i,j) = other(i,j);
			}
		}
		
		return *this;
	}


	size_t rows() const {
		return static_cast<const Derived*>(mTransposeSource)->cols();
	}
	size_t cols() const {
		return static_cast<const Derived*>(mTransposeSource)->rows();
	}

	const ScalarType& operator()(const size_t& i, const size_t& j) const {
		return static_cast<const Derived*>(mTransposeSource)->operator()(j, i);
	}
	ScalarType& operator()(const size_t& i, const size_t& j) {
		return static_cast<Derived*>(mTransposeSource)->operator()(j, i);
	}
};

template <typename Derived, typename ScalarType, int NumRows, int NumCols>
struct Block : public MatrixBase<Block<Derived, ScalarType, NumRows, NumCols>, ScalarType, NumRows, NumCols> {
	typedef Block<Derived, ScalarType, NumRows, NumCols> matrix_type;

	Derived* mBlockSource;
	int row_index;
	int col_index;

	Block(Derived* block_source, int row_index, int col_index) :
		mBlockSource(block_source),
		row_index(row_index),
		col_index(col_index)
	{ }

	Block(const Block &other) :
		mBlockSource(other.mBlockSource),
		row_index(other.row_index),
		col_index(other.col_index)
	{ }

	template <typename OtherDerived>
	Derived& operator=(const OtherDerived& other) {
		unsigned int i,j;
		for (i = 0; i < rows(); i++) {
			for (j = 0; j < cols(); j++) {
				this->operator()(i,j) = other(i,j);
			}
		}
		
		return *mBlockSource;
	}
	
	template <typename OtherDerived>
	Dynamic<ScalarType, -1, -1> operator*(const OtherDerived& other) {
		Dynamic<ScalarType, -1, -1> result (Dynamic<ScalarType, -1, -1>::Zero(rows(), other.cols()));

		unsigned int i,j,k;
		for (i = 0; i < rows(); i++) {
			for (j = 0; j < other.cols(); j++) {
				for (k = 0; k < other.rows(); k++) {
					result (i,j) += operator()(i,k) * other(k,j);
				}
			}
		}
		return result;
	}


	size_t rows() const {
		return NumRows;
	}
	size_t cols() const {
		return NumCols;
	}

	const ScalarType& operator()(const size_t& i, const size_t& j) const {
		return static_cast<const Derived*>(mBlockSource)->operator()(row_index + i, col_index + j);
	}
	ScalarType& operator()(const size_t& i, const size_t& j) {
		return static_cast<Derived*>(mBlockSource)->operator()(row_index + i,col_index + j);
	}

	template <typename OtherDerived>
	Dynamic<ScalarType, -1, -1> operator+(const OtherDerived& other) {
		Dynamic<ScalarType, -1, -1> result (*this);
		result += other;
		return result;
	}

	private:
	Block() { assert(0 && "Invalid call!"); };

	ScalarType* data() {
		assert("invalid call");
		return NULL;
	}

	const ScalarType* data() const {
		assert("invalid call");
		return NULL;
	}

};

template <typename val_type, int NumRows, int NumCols>
struct Fixed : public MatrixBase<Fixed<val_type, NumRows, NumCols>, val_type, NumRows, NumCols> {
	typedef Fixed<val_type, NumRows, NumCols> matrix_type;
	typedef val_type value_type;

	val_type mData[NumRows * NumCols];

	Fixed() {
		for (size_t i = 0; i < NumRows * NumCols; i++) {
			mData[i] = static_cast<val_type> (0.);
		}
	}

	Fixed (const Fixed &other) {
		assert (NumRows == other.rows() && NumCols == other.cols() && "Error: matrix dimensions do not match!");
		memcpy (mData, other.data(), sizeof (val_type) * rows() * cols());
	}

	template <typename OtherDerived>
	Fixed(const OtherDerived &other) {
		std::cout << "FCC" << std::endl;
		assert (other.rows() == NumRows && other.cols() == NumCols);

		for (size_t i = 0; i < rows(); i++) {
			for (size_t j = 0; j < cols(); j++) {
				this->operator()(i,j) = other(i,j);
			}
		}
	}

	template <typename OtherDerived>
	Fixed& operator=(const OtherDerived& other) {
		if (static_cast<const void*>(this) != static_cast<const void*>(&other)) {
			std::cout << "Fixed AO" << std::endl;

			for (size_t i = 0; i < other.rows(); i++) {
				for (size_t j = 0; j < other.cols(); j++) {
					this->operator()(i,j) = other(i,j);
				}
			}
		}

		return *this;
	}

	Fixed (int rows, int cols) {
		assert (rows == NumRows);
		assert (cols == cols);
	}

	template <typename OtherDerived>
	Fixed& operator+=(const OtherDerived& other) {
		assert (NumRows == other.rows() && NumCols == other.cols() && "Error: matrix dimensions do not match!");

		for (size_t i = 0; i < rows(); i++) {
			for (size_t j = 0; j < cols(); j++) {
				this->operator()(i,j) += other(i,j);
			}
		}
		return *this;
	}

	val_type& operator[](const size_t& i) {
		return mData[i];
	}

	val_type& operator()(const size_t& i, const size_t& j) {
		return mData[i*NumCols + j];
	}

	val_type* data() {
		return mData;
	}

	const val_type* data() const {
		return mData;
	}

	const val_type& operator()(const size_t& i, const size_t& j) const {
		return mData[i*NumCols + j];
	}

	size_t cols() const {
		return NumCols;
	}

	size_t rows() const {
		return NumRows;
	}
};

template <typename ScalarType, int Rows = -1, int Cols = -1>
struct Dynamic : public MatrixBase<Dynamic<ScalarType, -1, -1>, ScalarType, -1, -1> {
	typedef Dynamic<ScalarType> matrix_type;

	ScalarType* mData;
	int mNumRows;
	int mNumCols;

	Dynamic() : 
		mData (NULL),
		mNumRows (0),
		mNumCols (0)
	{
	}

	Dynamic(int num_rows, int num_cols) :
		mNumRows (num_rows),
		mNumCols (num_cols) {
		mData = new ScalarType[mNumRows * mNumCols];

		for (size_t i = 0; i < (size_t) mNumRows * (size_t) mNumCols; i++) {
			mData[i] = static_cast<ScalarType> (0.);
		}
	}

	~Dynamic() {
		if (mData != NULL) {
			delete[] mData;
			mData = NULL;
		}
	}

	Dynamic(const Dynamic &other) :
		mNumRows (other.rows()),
		mNumCols (other.cols()) {
		mData = new ScalarType[mNumRows * mNumCols];
		assert (mNumRows == (int) other.rows() && mNumCols == (int) other.cols() && "Error: matrix dimensions do not match!");
		memcpy (mData, other.data(), sizeof (ScalarType) * rows() * cols());
	}

	template <typename OtherDerived>
	Dynamic(const OtherDerived &other) :
		mNumRows (other.rows()),
		mNumCols (other.cols()) {
		mData = new ScalarType[mNumRows * mNumCols];
		assert (mNumRows == other.rows() && mNumCols == other.cols() && "Error: matrix dimensions do not match!");
		for (size_t i = 0; i < rows(); i++) {
			for (size_t j = 0; j < cols(); j++) {
				this->operator()(i,j) = other(i,j);
			}
		}
	}

	Dynamic& operator=(const Dynamic &other) {
		if (static_cast<const void*>(this) != static_cast<const void*>(&other)) {
			assert (mData != NULL);
			delete[] mData;

			mNumRows = other.rows();
			mNumCols = other.cols();
			mData = new ScalarType[mNumRows * mNumCols];
			
			for (size_t i = 0; i < rows(); i++) {
				for (size_t j = 0; j < cols(); j++) {
					this->operator()(i,j) = other(i,j);
				}
			}
		}

		return *this;
	}

	template <typename OtherDerived>
	Dynamic& operator=(const OtherDerived &other) {
		if (static_cast<const void*>(this) != static_cast<const void*>(&other)) {
			assert (mData != NULL);
			delete[] mData;

			mNumRows = other.rows();
			mNumCols = other.cols();
			mData = new ScalarType[mNumRows * mNumCols];
			
			for (size_t i = 0; i < rows(); i++) {
				for (size_t j = 0; j < cols(); j++) {
					this->operator()(i,j) = other(i,j);
				}
			}
		}

		return *this;
	}

	template <typename OtherDerived>
	Dynamic& operator+=(const OtherDerived& other) {
		assert (mNumRows == (int) other.rows() && mNumCols == (int) other.cols() && "Error: matrix dimensions do not match!");

		for (size_t i = 0; i < rows(); i++) {
			for (size_t j = 0; j < cols(); j++) {
				this->operator()(i,j) += other(i,j);
			}
		}
		
		return *this;
	}

	ScalarType& operator[](const size_t& i) {
		return mData[i];
	}

	ScalarType& operator()(const size_t& i, const size_t& j) {
		return mData[i*cols() + j];
	}

	const ScalarType& operator()(const size_t& i, const size_t& j) const {
		return mData[i*cols() + j];
	}

	ScalarType* data() {
		return mData;
	}

	const ScalarType* data() const {
		return mData;
	}

	size_t cols() const {
		return mNumCols;
	}

	size_t rows() const {
		return mNumRows;
	}
};

template <typename Derived, typename ScalarType, int Rows, int Cols>
inline std::ostream& operator<<(std::ostream& output, const MatrixBase<Derived, ScalarType, Rows, Cols> &matrix) {
	size_t max_width = 0;
	size_t out_width = output.width();

	// get the widest number
	for (size_t i = 0; i < matrix.rows(); i++) {
		for (size_t j = 0; j < matrix.cols(); j++) {
			std::stringstream out_stream;
			out_stream << matrix(i,j);
			max_width = std::max (out_stream.str().size(),max_width);
		}
	}

	// overwrite width if it was explicitly prescribed
	if (out_width != 0) {
		max_width = out_width;
	}

	for (unsigned int i = 0; i < matrix.rows(); i++) {
		output.width(0);
		output << "[ ";
		output.width(out_width);
		for (unsigned int j = 0; j < matrix.cols(); j++) {
			std::stringstream out_stream;
			out_stream.width (max_width);
			out_stream << matrix(i,j);
			output << out_stream.str();

			if (j < matrix.cols() - 1)
				output << ", ";
		}
		output << " ]";
		
		if (matrix.rows() > 1 && i < matrix.rows() - 1)
			output << std::endl;
	}
	return output;
}

}
