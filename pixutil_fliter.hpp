#if !defined(_BTK_PROJECT_PIXUTIL_FLITER_)
#define _BTK_PROJECT_PIXUTIL_FLITER_
#include "pixutil.hpp"
#include <cstdlib>
#include <cmath>

#ifndef PIXUTIL_MALLOC
    #define PIXUTIL_REALLOC std::realloc
    #define PIXUTIL_MALLOC std::malloc
    #define PIXUTIL_FREE std::free
#endif

namespace PixFliter{
namespace _Math{
    constexpr double PI = 3.14159265358979323846;
    constexpr double PI_2 = PI / 2;
    constexpr double PI_4 = PI / 4;
    constexpr double PI_8 = PI / 8;
} // namespace _Math
} // namespace PixFliter

namespace PixFliter{
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
} // namespace _Mem
} // namespace PixFliter


namespace PixFliter{
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
}

#endif // _BTK_PROJECT_PIXUTIL_FLITER_
