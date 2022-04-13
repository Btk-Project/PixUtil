#if !defined(_BTK_PROJECT_PIXUTIL_FILTER_)
#define _BTK_PROJECT_PIXUTIL_FILTER_
#include <initializer_list>
#include <complex>
#include <cstdlib>
#include <cmath>

#include "pixutil.hpp"

namespace PixFilter{
namespace _Math{
    constexpr double PI = 3.14159265358979323846;
    constexpr double PI_2 = PI / 2;
    constexpr double PI_4 = PI / 4;
    constexpr double PI_8 = PI / 8;
    constexpr double E = 2.71828182845904523536;
} // namespace _Math
} // namespace PixFilter

namespace PixFilter{
namespace _Mem{
    //< For memory allocation and deallocation
    struct AllocBase{
        void *operator new(size_t size){
            return PIXUTIL_MALLOC(size);
        }
        void operator delete(void *p){
            PIXUTIL_FREE(p);
        }

        void *operator new[](size_t size){
            return PIXUTIL_MALLOC(size);
        }
        void operator delete[](void *p){
            PIXUTIL_FREE(p);
        }
    };
    //< Internal simplest dynamic array for pod type
    template<typename T>
    class Vector:public AllocBase{
        public:
            using value_type = T;
        public:

            Vector(size_t size):
                _size(size),
                _data(static_cast<T*>(PIXUTIL_MALLOC(sizeof(T) * size)))
            {}
            Vector(size_t size,T v):
                _size(size),
                _data(static_cast<T*>(PIXUTIL_MALLOC(sizeof(T) * size)))
            {
                for(size_t i = 0;i < size;++i){
                    _data[i] = v;
                }
            }
            Vector(const Vector &) = delete;
            Vector(Vector &&vec):
                _size(vec._size),
                _data(vec._data){
                    
                vec._data = nullptr;
            }
            ~Vector(){
                PIXUTIL_FREE(_data);
            }
            T &operator[](size_t index){
                PIXUTIL_ASSERT(index < _size);
                return _data[index];
            }
            const T &operator[](size_t index) const{
                PIXUTIL_ASSERT(index < _size);
                return _data[index];
            }
            T &at(size_t index){
                PIXUTIL_ASSERT(index < _size);
                return _data[index];
            }
            const T &at(size_t index) const{
                PIXUTIL_ASSERT(index < _size);
                return _data[index];
            }
            size_t size() const{
                return _size;
            }
            void resize(size_t size){
                if(size == _size){
                    return;
                }
                T *new_data = static_cast<T*>(PIXUTIL_REALLOC(_data,sizeof(T) * size));
                if(new_data != nullptr){
                    _data = new_data;
                    _size = size;
                }
                else{
                    //TODO: throw exception or something
                }
            }

            T *data(){
                return _data;
            }
            const T *data() const{
                return _data;
            }

            T *begin(){
                return _data;
            }
            const T *begin() const{
                return _data;
            }
            T *end(){
                return _data + _size;
            }
            const T *end() const{
                return _data + _size;
            }
            void zero(){
                PIXUTIL_MEMSET(_data,0,sizeof(T) * _size);
            }
        private:
            size_t _size;
            T *_data;
        template<typename U>
        friend class Matrix;
    };
    //TODO add MatrixExpr to implement lazy evaluation
    template<typename T,size_t Row,size_t Col>
    class MatrixBlock{
        public:
            MatrixBlock(T *mat,size_t row,size_t col):
                mat(mat),
                cur_row(row),
                cur_col(col)
            {}
            MatrixBlock(const MatrixBlock &) = default;
            MatrixBlock(MatrixBlock &&) = default;
            ~MatrixBlock() = default;


            MatrixBlock &operator =(const T &v){
                //Assert v row and col
                PIXUTIL_ASSERT(v.row() == Row);
                PIXUTIL_ASSERT(v.col() == Col);
                for(size_t i = 0;i < Row;++i){
                    for(size_t j = 0;j < Col;++j){
                        mat->at(cur_row + i,cur_col + j) = v.at(i,j);
                    }
                }
                return *this;
            }
            template<class ...Args>
            MatrixBlock &assign(Args &&...args){
                static_assert(Row * Col == sizeof...(args),"num of Elem doesnot matched");
                auto l = {args...};
                
                size_t row = 0;
                size_t col = 0;

                for(auto v:l){
                    mat->at(row + cur_row,col + cur_col) = v;
                    col += 1;
                    if(col == Col){
                        col = 0;
                        row += 1;
                    }
                }
            }
        private:
            T *mat;
            size_t cur_row = 0;
            size_t cur_col = 0;
    };
    /**
     * @brief Simplest matrix
     * 
     * @tparam T 
     */
    template<typename T>
    class Matrix:public AllocBase{
        public:
            using value_type = T;
        public:
            /**
             * @brief Construct a new Matrix object
             * 
             */
            Matrix():
                _row(0),
                _col(0),
                _data(nullptr){
            }
            /**
             * @brief Construct a new Matrix object
             * 
             * @param row 
             * @param col 
             */
            Matrix(size_t row,size_t col):
                _row(row),
                _col(col),
                _data(static_cast<T*>(PIXUTIL_MALLOC(sizeof(T) * row * col)))
            {}
            /**
             * @brief Construct a new Matrix object from vector
             * 
             * @param row 
             * @param col 
             * @param vec 
             */
            Matrix(size_t row,size_t col,const Vector<T> &vec):
                _row(row),
                _col(col),
                _data(static_cast<T*>(PIXUTIL_MALLOC(sizeof(T) * row * col)))
            {
                PIXUTIL_ASSERT(row * col == vec.size());
                PIXUTIL_MEMCPY(_data,vec.data(),sizeof(T) * row * col);
            }
            /**
             * @brief Construct a new Matrix object from vector
             * 
             * @param row 
             * @param col 
             * @param vec 
             */
            Matrix(size_t row,size_t col,Vector<T> &&vec):
                _row(row),
                _col(col),
                _data(vec.data()){

                PIXUTIL_ASSERT(row * col == vec.size());
                vec._data = nullptr;
                vec._size = 0;
            }
            /**
             * @brief Construct a new Matrix object
             * 
             * @param mat 
             */
            Matrix(const Matrix & mat):Matrix(mat._row,mat._col){
                PIXUTIL_MEMCPY(_data,mat._data,sizeof(T) * _row * _col);
            }
            /**
             * @brief Construct a new Matrix object by different type
             * 
             * @tparam Elem 
             * @param mat 
             */
            template<typename Elem>
            Matrix(const Matrix<Elem> & mat):Matrix(mat._row,mat._col){
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < _col;++j){
                        _data[i * _col + j] = static_cast<T>(mat._data[i * _col + j]);
                    }
                }
            }

            Matrix(Matrix &&mat):
                _row(mat._row),
                _col(mat._col),
                _data(mat._data){
                    
                mat._data = nullptr;
            }
            ~Matrix(){
                PIXUTIL_FREE(_data);
            }
            T &at(size_t row,size_t col) noexcept{
                PIXUTIL_ASSERT(row < _row);
                PIXUTIL_ASSERT(col < _col);
                return _data[row * _col + col];
            }
            const T &at(size_t row,size_t col) const noexcept{
                PIXUTIL_ASSERT(row < _row);
                PIXUTIL_ASSERT(col < _col);
                return _data[row * _col + col];
            }
            size_t row() const noexcept{
                return _row;
            }
            //< For compatibility with View
            size_t height() const noexcept{
                return _row;
            }
            size_t col() const noexcept{
                return _col;
            }
            //< For compatibility with View
            size_t width() const noexcept{
                return _col;
            }
            /**
             * @brief Resize the Matrix
             * 
             * @param row 
             * @param col 
             */
            void resize(size_t row,size_t col){
                if(row == _row && col == _col){
                    return;
                }
                T *new_data = static_cast<T*>(PIXUTIL_MALLOC(sizeof(T) * row * col));
                if(new_data != nullptr){
                    _data = new_data;
                    _row = row;
                    _col = col;
                }
                else{
                    //TODO: throw exception or something
                }
            }
            void zero(){
                PIXUTIL_MEMSET(_data,0,sizeof(T) * _row * _col);
            }

            T *data(){
                return _data;
            }
            const T *data() const{
                return _data;
            }
            /**
             * @brief Compare two Matrix
             * 
             * @param mat 
             * @return true 
             * @return false 
             */
            bool compare(const Matrix &mat) const{
                if(_row != mat._row || _col != mat._col){
                    return false;
                }
                return PIXUTIL_MEMCMP(_data,mat._data,sizeof(T) * _row * _col) == 0;
            }
            template<typename Elem>
            bool compare(const Matrix<Elem> &mat) const{
                if(_row != mat._row || _col != mat._col){
                    return false;
                }
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < _col;++j){
                        if(_data[i * _col + j] != static_cast<T>(mat._data[i * _col + j])){
                            return false;
                        }
                    }
                }
                return true;
            }

            //Operator
            T *operator [](size_t row){
                PIXUTIL_ASSERT(row < _row);
                return _data + row * _col;
            }
            const T *operator [](size_t row) const{
                PIXUTIL_ASSERT(row < _row);
                return _data + row * _col;
            }
            //Assignment
            Matrix &operator =(const Matrix &mat){
                if(this == &mat){
                    return *this;
                }
                if(_row != mat._row || _col != mat._col){
                    resize(mat._row,mat._col);
                }
                PIXUTIL_MEMCPY(_data,mat._data,sizeof(T) * _row * _col);
                return *this;
            }
            Matrix &operator =(Matrix &&mat){
                if(this == &mat){
                    return *this;
                }
                PIXUTIL_FREE(_data);
                _row = mat._row;
                _col = mat._col;
                _data = mat._data;
                mat._data = nullptr;
                mat._row = 0;
                mat._col = 0;
                return *this;
            }
        public:
            //Math operators

            //For number
            template<typename Elem>
            Matrix &operator +=(Elem num){
                for(size_t i = 0;i < _row * _col;++i){
                    _data[i] += num;
                }
                return *this;
            }
            template<class Elem>
            Matrix &operator -=(Elem num){
                for(size_t i = 0;i < _row * _col;++i){
                    _data[i] -= num;
                }
                return *this;
            }
            template<class Elem>
            Matrix &operator *=(Elem num){
                for(size_t i = 0;i < _row * _col;++i){
                    _data[i] *= num;
                }
                return *this;
            }
            template<class Elem>
            Matrix &operator /=(Elem num){
                for(size_t i = 0;i < _row * _col;++i){
                    _data[i] /= num;
                }
                return *this;
            }
            //For Matrix
            Matrix &operator +=(const Matrix &mat){
                PIXUTIL_ASSERT(_row == mat._row && _col == mat._col);
                for(size_t i = 0;i < _row * _col;++i){
                    _data[i] += mat._data[i];
                }
                return *this;
            }
            Matrix &operator -=(const Matrix &mat){
                PIXUTIL_ASSERT(_row == mat._row && _col == mat._col);
                for(size_t i = 0;i < _row * _col;++i){
                    _data[i] -= mat._data[i];
                }
                return *this;
            }

            template<typename Elem>
            Matrix operator +(Elem num) const{
                Matrix<T> mat(*this);
                mat += num;
                return mat;
            }
            template<typename Elem>
            Matrix operator -(Elem num) const{
                Matrix<T> mat(*this);
                mat -= num;
                return mat;
            }
            template<typename Elem>
            Matrix operator *(Elem num) const{
                Matrix<T> mat(*this);
                mat *= num;
                return mat;
            }
            template<typename Elem>
            Matrix operator /(Elem num) const{
                Matrix<T> mat(*this);
                mat /= num;
                return mat;
            }
            //Dot
            Matrix &operator *=(const Matrix &mat){
                PIXUTIL_ASSERT(_col == mat._row);
                //OMP Optimize
                Matrix<T> tmp(_row,mat._col);
                PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(tmp)")
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < mat._col;++j){
                        tmp.at(i,j) = 0;
                        for(size_t k = 0;k < _col;++k){
                            tmp.at(i,j) += at(i,k) * mat.at(k,j);
                        }
                    }
                }
                *this = std::move(tmp);
                return *this;
            }
            Matrix operator *(const Matrix &mat) const{
                PIXUTIL_ASSERT(_col == mat._row);
                //OMP Optimize
                Matrix<T> tmp(_row,mat._col);
                PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(tmp)")
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < mat._col;++j){
                        tmp.at(i,j) = 0;
                        for(size_t k = 0;k < _col;++k){
                            tmp.at(i,j) += at(i,k) * mat.at(k,j);
                        }
                    }
                }
                return tmp;
            }

            //Transpose
            Matrix transpose() const{
                Matrix<T> tmp(_col,_row);
                PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(tmp)")
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < _col;++j){
                        tmp.at(j,i) = _data[i * _col + j];
                    }
                }
                return tmp;
            }
            //Inverse
            Matrix inverse() const{
                PIXUTIL_ASSERT(_row == _col);
                Matrix<T> tmp(_row,_col);
                PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(tmp)")
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < _col;++j){
                        tmp.at(i,j) = _data[j * _row + i];
                    }
                }
                return tmp;
            }
            /**
             * @brief Map the matrix to a new matrix
             * 
             * @tparam Fn 
             * @param fn 
             * @return Matrix 
             */
            template<typename Fn>
            Matrix map(Fn &&fn) const{
                Matrix<T> tmp(_row,_col);
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < _col;++j){
                        tmp.at(i,j) = fn(_data[i * _col + j]);
                    }
                }
                return tmp;
            }
            /**
             * @brief Apply a function to each element of the matrix
             * 
             * @tparam Fn 
             * @param fn 
             * @return Matrix& 
             */
            template<typename Fn>
            Matrix &apply(Fn &&fn){
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < _col;++j){
                        _data[i * _col + j] = fn(_data[i * _col + j]);
                    }
                }
                return *this;
            }
        public:
            //Helper / cast
            template<size_t Row,size_t Col>
            MatrixBlock<Matrix,Row,Col> block(size_t r,size_t c){
                PIXUTIL_ASSERT(r + Row <= _row && c + Col <= _col);
                return MatrixBlock<Matrix,Row,Col>(this,r,c);
            }


            template<typename Elem>
            Matrix<Elem> cast() const{
                return *this;
            }

            operator Vector<T>() const{
                Vector<T> vec(_row * _col);
                PIXUTIL_MEMCPY(vec.data(),_data,sizeof(T) * _row * _col);
                return vec;
            }

            template<typename Elem>
            operator Vector<Elem>() const{
                Vector<Elem> vec(_row * _col);

                for(size_t i = 0;i < _row * _col;++i){
                    vec.at(i) = static_cast<Elem>(_data[i]);
                }
                
                return vec;
            }
        private:
            size_t _row;
            size_t _col;
            T *_data;
    };
    template<typename T,size_t Row,size_t Col>
    class FixedMatrix:public Matrix<T>{
        public:
            FixedMatrix():Matrix<T>(Row,Col){};
    };

    #ifdef PIXUTIL_MATRIX_IO
    /**
     * @brief Print matrix
     * 
     * @tparam T 
     * @param os 
     * @param mat 
     * @return std::ostream& 
     */
    template<typename T>
    std::ostream &operator <<(std::ostream &os,const Matrix<T> &mat){
        for(size_t i = 0;i < mat.row();++i){
            os << "[";
            for(size_t j = 0;j < mat.col();++j){
                os << mat.at(i,j) << " ";
                if(j != mat.col() - 1){
                    os << ", ";
                }
            }
            os << "]" << std::endl;
        }
        return os;
    }
    #endif
} // namespace _Mem
} // namespace PixFilter


namespace PixFilter{
    using TransformMatrix = _Mem::Matrix<double>;
    template<typename T>
    using Matrix = _Mem::Matrix<T>;
    template<typename T>
    using Vector = _Mem::Vector<T>;
    template<typename ...Args>
    using View = PixUtil::View<Args...>;
    using Color = PixUtil::Color;
    //Useful typedef
    using Uint16 = PixUtil::Uint16;
    using Uint32 = PixUtil::Uint32;
    using Uint64 = PixUtil::Uint64;
    using Uint8 = PixUtil::Uint8;
    //Useful function
    using PixUtil::clamp;
    using PixUtil::min;
    using PixUtil::max;

    //TODO: GaussianBlur
    namespace _Math{
        template<typename T>
        inline bool IsPowerOfTwo(T n){
            //Kepp dividing 2
            return n && !(n & (n - 1));
        }
        template<typename T>
        inline T    ToPowerOfTwo(T n){
            //Round up to the nearest power of 2
            T power = 1;
            while(power < n){
                power *= 2;
            }
            return power;
        }
        inline void GaussianFliter(Vector<float> &fliter,float sigma,int r){
            float sum = 0;
            for (int i = 0;i<r;i++){
                for(int j = 0;j <r;j++){
                    float r2 = (i - r/2)*(i - r/2) + (j - r/2)*(j - r/2);
                    fliter.at(i*r + j) = (1.0 / (2* PI * sigma * sigma)) * std::exp(-r2 / (2* sigma * sigma));
                    sum += fliter.at(i*r + j);
                }
            }
            int i = 0;
            switch ((r*r)%8){
                for(;i<r*r;){
                    case 0:
                        fliter.at(i) /= sum;i++;
                    case 7:
                        fliter.at(i) /= sum;i++;
                    case 6:
                        fliter.at(i) /= sum;i++;
                    case 5:
                        fliter.at(i) /= sum;i++;
                    case 4:
                        fliter.at(i) /= sum;i++;
                    case 3:
                        fliter.at(i) /= sum;i++;
                    case 2:
                        fliter.at(i) /= sum;i++;
                    case 1:
                        fliter.at(i) /= sum;i++;
                }
            }
        }
        //2D Fast Fourier Transform type
        using FFTMatrix = Matrix<std::complex<double>>;
        using FFTVector = Vector<std::complex<double>>;
        using FFTComplex = std::complex<double>;
        /**
         * @brief Recursive 1D FFT function(from github wareya/fft)
         * 
         * @param outputs 
         * @param inputs 
         * @param size 
         * @param gap must be 1
         * @param inverse 
         */
        inline void FFT1DCore(FFTComplex *outputs,const FFTComplex *inputs,size_t size,size_t gap,bool inverse){
            #if 0
            //Use the recursive version
            if(size == 1){
                outputs[0] = inputs[0];
                return;
            }
            FFT1DCore(outputs           ,inputs      ,size / 2,gap * 2,inverse);
            FFT1DCore(outputs + size / 2,inputs + gap,size / 2,gap * 2,inverse);
            for(size_t i = 0;i < size / 2;++i){
                //Begin with the twiddle factor
                FFTComplex t = outputs[i           ];
                FFTComplex u = outputs[i + size / 2];

                //Twiddle factor
                FFTComplex twiddle = {
                    std::cos(2 * PI * i / size),
                    std::sin(2 * PI * i / size) * (inverse ? 1 : -1)
                };
                double bias_real = u.real() * twiddle.real() - u.imag() * twiddle.imag();
                double bias_imag = u.imag() * twiddle.real() + u.real() * twiddle.imag();

                //Output
                outputs[i           ] = t + FFTComplex{bias_real,bias_imag};
                outputs[i + size / 2] = t - FFTComplex{bias_real,bias_imag};
            }
            #else
            //Use the iterative version
            PIXUTIL_MEMCPY(outputs,inputs,sizeof(FFTComplex) * size);
            PIXUTIL_UNUSED(gap);

            int lim = int(size);
            int len = 0;
            while((1 << len) < lim){
                ++ len;
            }

            Vector<int> bit(lim,0);
            for(int i = 0;i < lim;++i){
                bit.at(i) = (bit.at(i >> 1) >> 1) | ((i & 1) << (len - 1));
            }

            for(int i = 0;i < lim;++i){
                if(i < bit.at(i)){
                    FFTComplex t = outputs[i];
                    outputs[i] = outputs[bit.at(i)];
                    outputs[bit.at(i)] = t;
                }
            }

            int opt = (inverse) ? 1 : -1;
            for(int m = 1;m <= lim;m <<= 1){
                int mh = m >> 1;
                FFTComplex wm = {
                    std::cos(2 * PI / m),
                    std::sin(2 * PI / m) * opt
                };
                for(int i = 0;i < lim;i += m){
                    FFTComplex w = {1,0};
                    for(int j = i;j < i + mh;++j){
                        FFTComplex t = outputs[j + mh] * w;
                        FFTComplex u = outputs[j];
                        outputs[j] = u + t;
                        outputs[j + mh] = u - t;
                        w *= wm;
                    }
                }
            }
            #endif
        }
        inline void FFT1D(FFTComplex *outputs,const FFTComplex *inputs,size_t size,bool inverse){
            FFT1DCore(outputs,inputs,size,1,inverse);
        }
        //2D
        inline void FFT2D(FFTMatrix &output,const FFTMatrix &input,bool inverse = false){
            //2D FFT to 1D FFT
            //First. Do 1D FFT on each row
            //Second. Do 1D FFT on each column

            //Check input row and col is power of 2
            size_t row = input.row();
            size_t col = input.col();

            PIXUTIL_ASSERT(IsPowerOfTwo(row));
            PIXUTIL_ASSERT(IsPowerOfTwo(col));
            
            //Resize output
            output.resize(row,col);
            //Do the first dimension
            for(size_t i = 0;i < row;++i){
                FFT1D(output[i],input[i],col,inverse);
            }
            //Do the second dimension

            //TODO: Optimize here,avoid copy
            FFTVector input_tmp(row);
            FFTVector output_tmp(row);
            for(size_t i = 0;i < col;++i){
                for(size_t j = 0;j < row;++j){
                    input_tmp[j] = output[j][i];
                }
                FFT1D(output_tmp.data(),input_tmp.data(),row,inverse);
                for(size_t j = 0;j < row;++j){
                    output[j][i] = output_tmp[j];
                }
            }
            //Half the result
            output /= std::sqrt(row * col);
        }
        inline void IFFT2D(FFTMatrix &output,const FFTMatrix &input){
            FFT2D(output,input,true);
        }
        /**
         * @brief Shift 
         * 
         * @param mat 
         */
        inline void FFTShift(FFTMatrix &mat){
            //Shift to center
            for(size_t i = 0;i < mat.row();++i){
                for(size_t j = 0;j < mat.col();++j){
                    mat.at(i,j) *= std::pow(-1,i + j);
                }
            }
        }
        inline void DFT(FFTMatrix &dst,const FFTMatrix &src,bool inverse = false){
            size_t w = src.col();
            size_t h = src.row();

            FFTMatrix mat_w(w,w);
            FFTMatrix mat_h(h,h);
            //Check is inversed
            double v;
            if(inverse){
                v = 2;
            }
            else{
                v = -2;
            }
            //Prepare transform matrix
            PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(mat_w,mat_h)")
            for(size_t i = 0;i < w;++i){
                for(size_t j = 0;j < w;++j){
                    mat_w.at(i,j) = std::exp(std::complex<double>(0,v * PI * i * j / w));
                }
            }
            PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(mat_w,mat_h)")
            for(size_t i = 0;i < h;++i){
                for(size_t j = 0;j < h;++j){
                    mat_h.at(i,j) = std::exp(std::complex<double>(0,v * PI * i * j / h));
                }
            }
            //Transform
            dst  = mat_h * src * mat_w;
            dst /= std::sqrt(w * h);
        }
        inline void IDFT(FFTMatrix &dst,const FFTMatrix &src){
            DFT(dst,src,true);
        }
        inline Uint8 ComplexToUint8(const std::complex<double> &c) noexcept{
            return clamp(std::abs(c),0.0,255.0);
        }
    } // namespace _Math
    using _Math::FFTMatrix;
    using _Math::FFTShift;
    using _Math::IDFT;
    using _Math::DFT;

    template<typename View>
    void GaussianBlur(View &view,float sigma,int radius){
        if(sigma < 3/2 && radius < 1){
            return;
        }

        int w = view.width();
        int h = view.height();

        //Alloc a temporary buffer
        _Mem::Vector<uint8_t> image(
            w * h * view.bytes_per_pixel()
        );
        View buf_view(image.data(),w,h,view.traits());
        //Alloc fliter
        int blur_size = static_cast<int>(sigma * 6) + 1;
        if(radius > 0) blur_size = radius;
        //Build fliter
        _Mem::Vector<float> fliter(blur_size * blur_size);
        _Math::GaussianFliter(fliter,sigma,radius);

        //Blur
        PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(buf_view,fliter)")
        for(int i = 0;i < h;i ++){
            for(int j = 0;j < w;j ++){
                float sum[4] = {0.0f,0.0f,0.0f,0.0f};
                int index = -1;

                for(int m = i - blur_size / 2;m < i + blur_size / 2; m ++){
                    for(int n = j - blur_size / 2;n < j + blur_size / 2;n ++) {
                        ++ index;

                        int x = n;
                        int y = m;

                        if(x < 0){
                            x = - x;
                        }
                        if(y < 0){
                            y = - y;
                        }
                        //Flip if x or y is out of range
                        if(x >= w){
                            x = 2 * w - x - 1;
                        }
                        if(y >= h){
                            y = 2 * h - y - 1;
                        }
                        Color color = view[y][x];

                        sum[0] += color.r * fliter.at(index);
                        sum[1] += color.g * fliter.at(index);
                        sum[2] += color.b * fliter.at(index);
                        sum[3] += color.a * fliter.at(index);
                    }
                }

                for(auto &v:sum){
                    v = clamp(v,0.0f,255.0f);
                }
                buf_view[i][j] = Color{
                    static_cast<Uint8>(sum[0]),
                    static_cast<Uint8>(sum[1]),
                    static_cast<Uint8>(sum[2]),
                    static_cast<Uint8>(sum[3])
                };
            }
        }

        //Copy back
        PixUtil::CopyPixels(view,buf_view);
    }
    /**
     * @brief Filter 
     * 
     * @tparam View 
     * @tparam T 
     * @param dst Output view
     * @param src Input view
     * @param filter Input matrix
     */
    template<typename View,typename T>
    void Filter2D(View &dst,const View &src,const Matrix<T> &filter){
        #ifndef PIXUTIL_NO_RUNTIME_CHECK
        if(&dst == &src){
            int w = src.width();
            int h = src.height();
            //Alloc a temporary buffer
            _Mem::Vector<uint8_t> image(
                w * h * src.bytes_per_pixel()
            );
            View buf_view(image.data(),w,h,src.traits());
            //Do filter again
            Filter2D(buf_view,src,filter);
            //Copy back
            PixUtil::CopyPixels(dst,buf_view);
            return;
        }
        #endif
        //TODO: refactor
        int height = src.height();
        int width = src.width();
        int filter_height = filter.row();
        int filter_width = filter.col();
        //Process one channal 
        PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(dst,src,filter)")
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                double sum[3] = {0.0,0.0,0.0};

                for (int m = 0; m < filter_height; m++) {
                    for (int n = 0; n < filter_width; n++) {
                        int x = i + m - filter_height / 2;
                        int y = j + n - filter_width / 2;

                        if (x < 0) {
                            x = -x;
                        }
                        if (y < 0) {
                            y = -y;
                        }
                        if (x >= height) {
                            x = 2 * height - x - 1;
                        }
                        if (y >= width) {
                            y = 2 * width - y - 1;
                        }

                        Color color = src[x][y];

                        sum[0] += color.r * filter.at(m, n);
                        sum[1] += color.g * filter.at(m, n);
                        sum[2] += color.b * filter.at(m, n);
                        // sum[3] += color.a * filter.at(m, n);
                    }
                }
                //FIXME:It will cause output whole black if we process the alpha channal
                for(auto &v:sum){
                    v = clamp(v,0.0,255.0);
                }
                dst[i][j] = Color{
                    static_cast<Uint8>(sum[0]),
                    static_cast<Uint8>(sum[1]),
                    static_cast<Uint8>(sum[2]),
                    src[i][j].to_color().a//So we use the original alpha
                };
            }
        }
    }
    //FFT 
    template<typename View>
    void MergeFFTMatrixToView(View &dst,const FFTMatrix inputs[],int n){
        PIXUTIL_ASSERT(n > 0 && n <= 4);
        #ifndef PIXUTIL_NO_RUNTIME_CHECK
        //Check the size
        if(inputs[0].row() != dst.height() || inputs[0].col() != dst.width()){
            //Create a temporary buffer and begin scale
            //TODO: Should we add a helper function to handle this?
            //TODO: May be fill matrix extra space in 0 is better?
            Vector<Uint8> image(
                inputs[0].row() * inputs[0].col() * dst.bytes_per_pixel()
            );
            View buf_view(image.data(),inputs[0].col(),inputs[0].row(),dst.traits());
            PixUtil::MergeChannels(buf_view,inputs,n,_Math::ComplexToUint8);
            //Scale back
            PixUtil::BilinearScale(dst,buf_view);
            return;
        }
        #endif
        PixUtil::MergeChannels(dst,inputs,n,_Math::ComplexToUint8);
    }
    /**
     * @brief FFT2D (It requires w and h is power of 2)
     * 
     * @tparam View 
     * @param outputs 
     * @param n 
     * @param src 
     */
    template<typename View>
    void FFT2D(_Math::FFTMatrix outputs[],int n,const View &src){
        #ifndef PIXUTIL_NO_RUNTIME_CHECK
        //Dymamic check
        size_t w = src.width();
        size_t h = src.height();
        if(!_Math::IsPowerOfTwo(w) || !_Math::IsPowerOfTwo(h)){
            //TODO: Should we add a helper function to handle this?
            w = _Math::ToPowerOfTwo(w);
            h = _Math::ToPowerOfTwo(h);

            Vector<Uint8> image(w * h * src.bytes_per_pixel());
            View buf_view(image.data(),w,h,src.traits());
            //Scale to power of 2
            PixUtil::BilinearScale(buf_view,src);
            //Calculate FFT
            return FFT2D(outputs,n,buf_view);
        }
        #else
        PIXUTIL_ASSERT(
            _Math::IsPowerOfTwo(src.width()) && _Math::IsPowerOfTwo(src.height())
        );
        #endif
        PIXUTIL_ASSERT(n > 0 && n <= 4);

        _Math::FFTMatrix *inputs = new _Math::FFTMatrix[n];
        //Pack the view to n matrix
        for(int i = 0;i < n;i++){
            inputs[i].resize(src.height(),src.width());
        }
        PixUtil::SplitChannels(inputs,n,src);
        //Begin FFT
        PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(inputs,outputs)")
        for(int i = 0;i < n;i ++){
            _Math::FFTShift(inputs[i]);
            _Math::FFT2D(outputs[i],inputs[i]);
        }
        delete [] inputs;
    }
    template<typename View,size_t N>
    void FFT2D(_Math::FFTMatrix (&outputs)[N],const View &src){
        static_assert(N > 0 && N <= 4,"N must be in [1,4]");
        FFT2D(outputs,N,src);
    }
    /**
     * @brief FFT2D a view to a view
     * 
     * @tparam View 
     * @param dst 
     * @param src 
     */
    template<typename View>
    void FFT2D(View &dst,const View &src){
        //We only process R G B channel
        _Math::FFTMatrix outputs[3];
        FFT2D(outputs,src);
        //Done in one step
        //Merge
        MergeFFTMatrixToView(dst,outputs,3);
    }
    /**
     * @brief Inverse FFT2D
     * 
     * @tparam View 
     * @param dst 
     * @param mats 
     * @param n 
     */
    template<typename View>
    void IFFT2D(View &dst,const _Math::FFTMatrix src[],int n){
        PIXUTIL_ASSERT(n > 0 && n <= 4);
        _Math::FFTMatrix *outputs = new _Math::FFTMatrix[n];
        //Init the output
        for(int i = 0;i < n;i++){
            outputs[i].resize(src[i].row(),src[i].col());
        }
        //Begin to process
        PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(src,outputs)")
        for(int i = 0;i < n;i++){
            _Math::IFFT2D(outputs[i],src[i]);
            _Math::FFTShift(outputs[i]);
        }
        //Map to dst
        MergeFFTMatrixToView(dst,outputs,n);
        delete [] outputs;
    }
    template<typename View,size_t N>
    void IFFT2D(View &dst,const _Math::FFTMatrix (&src)[N]){
        static_assert(N > 0 && N <= 4,"N must be in [1,4]");
        IFFT2D(dst,src,N);
    }

    /**
     * @brief DFT a view to three matrix
     * 
     * @tparam View 
     * @param outputs FFTMatrix[n]
     * @param n Number of matrix
     * @param src View
     */
    template<typename View>
    void DFT(_Math::FFTMatrix outputs[],int n,const View &src){
        PIXUTIL_ASSERT(n > 0 && n <= 4);

        _Math::FFTMatrix *inputs = new _Math::FFTMatrix[n];
        //Pack the view to three matrix
        for(int i = 0;i < n;i++){
            inputs[i].resize(src.height(),src.width());
        }
        PixUtil::SplitChannels(inputs,n,src);
        //Begin DFT
        PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(inputs,outputs)")
        for(int i = 0;i < n;i ++){
            _Math::FFTShift(inputs[i]);
            _Math::DFT(outputs[i],inputs[i]);
        }
        delete [] inputs;
    }
    /**
     * @brief Helper for DFT
     * 
     * @tparam View 
     * @tparam N 
     * @param src 
     */
    template<typename View,size_t N>
    void DFT(_Math::FFTMatrix (&outputs)[N],const View &src){
        static_assert(N > 0 && N <= 4,"N must be in [1,4]");
        DFT(outputs,N,src);
    }
    /**
     * @brief DFT a view to a view
     * 
     * @tparam View 
     * @param dst 
     * @param src 
     */
    template<typename View>
    void DFT(View &dst,const View &src){
        //We only process R G B channel
        _Math::FFTMatrix outputs[3];
        DFT(outputs,src);
        //Done in one step
        //Merge
        PixUtil::MergeChannels(dst,outputs,3,_Math::ComplexToUint8);
    }
    /**
     * @brief IDFT some matrix to a view
     * 
     * @tparam View 
     * @param dst 
     * @param src 
     * @param n 
     */
    template<typename View>
    void IDFT(View &dst,const _Math::FFTMatrix *src,int n){
        PIXUTIL_ASSERT(n > 0 && n <= 4);
        _Math::FFTMatrix *outputs = new _Math::FFTMatrix[n];
        //Init the output
        for(int i = 0;i < n;i++){
            outputs[i].resize(src[i].row(),src[i].col());
        }
        //Begin to process
        PIXUTIL_OMP_DECL("omp parallel for schedule(dynamic) shared(outputs,src)")
        for(int i = 0;i < n;i++){
            _Math::IDFT(outputs[i],src[i]);
            _Math::FFTShift(outputs[i]);
        }
        //Map to dst
        PixUtil::MergeChannels(dst,outputs,n,_Math::ComplexToUint8);
        delete [] outputs;
    }
    /**
     * @brief IDFT
     * 
     * @tparam View 
     * @param dst The output view(GrayView is a good choice)
     * @param mat 
     */
    template<typename View,size_t N>
    void IDFT(View &dst,const _Math::FFTMatrix (&mat)[N]){
        static_assert(N > 0 && N <= 4,"N must be in [1,4]");
        IDFT(dst,&mat,N);
    }
    //View Transform
    template<typename View1,typename View2>
    void Transform(View1 &dst,const View2 &src,const TransformMatrix &m){
        #if 0
        //Make src to dst mat
        auto mat = m.inverse();
        //Build done
        Matrix<double> cur_p(3,1);
        Matrix<double> src_p(3,1);

        for(int y = 0;y < dst.height();++ y){
            for(int x = 0;x < dst.width();++ x){
                cur_p.at(0,0) = y - center_y;
                cur_p.at(1,0) = x - center_x;
                cur_p.at(2,0) = 1;

                src_p = mat * cur_p;

                double src_y = src_p.at(0,0);
                double src_x = src_p.at(1,0);

                if(src.has_point(src_x,src_y)){
                    dst[y][x] = src[src_y][src_x].to_color();
                }
            }
        }
        #else
        PIXUTIL_UNUSED(m);
        PIXUTIL_UNUSED(src);
        PIXUTIL_UNUSED(dst);
        PIXUTIL_ASSERT(!"Not implemented");
        #endif
    }

    //Some useful function template

    /**
     * @brief Generate a matrix from a list
     * 
     * @tparam T 
     * @tparam Row 
     * @tparam Col 
     * @tparam Args 
     * @param args 
     * @return Matrix<T> 
     */
    template<typename T,size_t Row,size_t Col,typename ...Args>
    Matrix<T> ImportMatrix(Args &&...args){
        static_assert(Row * Col == sizeof...(args),"num of Elem doesnot matched");
        Matrix<T> mat(Row,Col);

        auto l = {std::forward<T>(args)...};

        size_t row = 0;
        size_t col = 0;

        for(auto v:l){
            mat.at(row,col) = v;
            col += 1;
            if(col == Col){
                col = 0;
                row += 1;
            }
        }

        return mat;
    }

    //some useful 2D Filter matrix 
    
    /**
     * @brief Sharpen filter
     * 
     */
    inline Matrix<double> SharpenFilter(){
        return ImportMatrix<double,3,3>(
            -1.0f, -1.0f, -1.0f,
            -1.0f,  9.0f, -1.0f,
            -1.0f, -1.0f, -1.0f
        );
    }
    /**
     * @brief Edge detect filter
     * 
     */
    inline Matrix<double> EdgeDetectFilter(){
        return ImportMatrix<double,3,3>(
            -1.0f, -1.0f, -1.0f,
            -1.0f,  8.0f, -1.0f,
            -1.0f, -1.0f, -1.0f
        );
    }
    /**
     * @brief Box blur filter
     * 
     * @param radius The radius of the box blur
     */
    inline Matrix<double> BoxBlurFilter(int radius){
        Matrix<double> mat(radius * 2 + 1,radius * 2 + 1);
        double sum = 0.0;
        for(size_t i = 0;i < mat.row();i++){
            for(size_t j = 0;j < mat.col();j++){
                mat.at(i,j) = 1.0 / (mat.row() * mat.col());
                sum += mat.at(i,j);
            }
        }
        return mat;
    }
}

#endif // _BTK_PROJECT_PIXUTIL_FILTER_
