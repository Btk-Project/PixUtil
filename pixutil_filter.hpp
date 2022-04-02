#if !defined(_BTK_PROJECT_PIXUTIL_FILTER_)
#define _BTK_PROJECT_PIXUTIL_FILTER_
#include <initializer_list>
#include <cstdlib>
#include <cmath>

#include "pixutil.hpp"

#ifndef PIXUTIL_MALLOC
    #define PIXUTIL_REALLOC std::realloc
    #define PIXUTIL_MALLOC std::malloc
    #define PIXUTIL_FREE std::free
#endif

#ifndef PIXUTIL_MEMCPY
    #define PIXUTIL_MEMCPY ::memcpy
#endif

namespace PixFilter{
namespace _Math{
    constexpr double PI = 3.14159265358979323846;
    constexpr double PI_2 = PI / 2;
    constexpr double PI_4 = PI / 4;
    constexpr double PI_8 = PI / 8;
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
        private:
            size_t _size;
            T *_data;
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
            T &at(size_t row,size_t col){
                PIXUTIL_ASSERT(row < _row);
                PIXUTIL_ASSERT(col < _col);
                return _data[row * _col + col];
            }
            const T &at(size_t row,size_t col) const{
                PIXUTIL_ASSERT(row < _row);
                PIXUTIL_ASSERT(col < _col);
                return _data[row * _col + col];
            }
            size_t row() const{
                return _row;
            }
            size_t col() const{
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

            T *data(){
                return _data;
            }
            const T *data() const{
                return _data;
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
            //Dot
            Matrix &operator *=(const Matrix &mat){
                PIXUTIL_ASSERT(_col == mat._row);
                Matrix<T> tmp(_row,mat._col);
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < mat._col;++j){
                        tmp.at(i,j) = 0;
                        for(size_t k = 0;k < _col;++k){
                            tmp.at(i,j) += _data[i * _col + k] * mat.at(k,j);
                        }
                    }
                }
                *this = std::move(tmp);
                return *this;
            }
            Matrix operator *(const Matrix &mat) const{
                PIXUTIL_ASSERT(_col == mat._row);
                Matrix<T> tmp(_row,mat._col);
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < mat._col;++j){
                        tmp.at(i,j) = 0;
                        for(size_t k = 0;k < _col;++k){
                            tmp.at(i,j) += _data[i * _col + k] * mat.at(k,j);
                        }
                    }
                }
                return tmp;
            }

            //Transpose
            Matrix transpose() const{
                Matrix<T> tmp(_col,_row);
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
                for(size_t i = 0;i < _row;++i){
                    for(size_t j = 0;j < _col;++j){
                        tmp.at(i,j) = _data[j * _row + i];
                    }
                }
                return tmp;
            }
        private:
            size_t _row;
            size_t _col;
            T *_data;
    };
    // template<typename T>
    // struct ImportMatrix{
    //     ImportMatrix(size_t row,size_t col):
    //         _row(row),
    //         _col(col){

    //     }

    //     size_t _row;
    //     size_t _col;
    // };

    // template<typename T,typename T2>
    // Matrix<T> operator <<(ImportMatrix<T> opt,std::initializer_list<T2> l){
    //     Matrix<T> mat(opt._row,opt._col);

    //     size_t row = 0;
    //     size_t col = 0;

    //     for(auto v:l){
    //         mat.at(row,col) = v;
    //         row += 1;
    //         if(row == opt._row){
    //             row = 0;
    //             col += 1;
    //         }
    //     }

    //     return mat;
    // }
} // namespace _Mem
} // namespace PixFilter


namespace PixFilter{
    template<typename T>
    using Matrix = _Mem::Matrix<T>;
    template<typename T>
    using Vector = _Mem::Vector<T>;
    template<typename ...Args>
    using View = PixUtil::View<Args...>;
    using Color = PixUtil::Color;

    //TODO: GaussianBlur
    namespace _Math{
        inline void GaussianFliter(_Mem::Vector<float> &fliter,float sigma,int r){
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
    }

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
                    v = PixUtil::clamp(v,0.0f,255.0f);
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
        int height = src.width();
        int width = src.height();
        int filter_height = filter.row();
        int filter_width = filter.col();
        //Process one channal 
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

                        Color color = src.at(x,y);

                        sum[0] += color.r * filter.at(m, n);
                        sum[1] += color.g * filter.at(m, n);
                        sum[2] += color.b * filter.at(m, n);
                        // sum[3] += color.a * filter.at(m, n);
                    }
                }
                //FIXME:It will cause output whole black if we process the alpha channal
                for(auto &v:sum){
                    v = PixUtil::clamp(v,0.0,255.0);
                }
                dst.at(i, j) = Color{
                    static_cast<Uint8>(sum[0]),
                    static_cast<Uint8>(sum[1]),
                    static_cast<Uint8>(sum[2]),
                    src.at(i, j).to_color().a//So we use the original alpha
                };
            }
        }
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
}

#endif // _BTK_PROJECT_PIXUTIL_FILTER_
