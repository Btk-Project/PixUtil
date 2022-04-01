#if !defined(_BTK_PROJECT_PIXUTIL_FLITER_)
#define _BTK_PROJECT_PIXUTIL_FLITER_
#include "pixutil.hpp"
#include <cmath>

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
        #ifdef PIXUTIL_MALLOC
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
        #endif
    };
    //< Internal simplest dynamic array for pod type
    template<typename T>
    class Vector:public AllocBase{
        public:
            Vector(size_t size):
                _size(size),
                _data(new T[size])
            {}
            ~Vector(){
                delete[] _data;
            }
            T &operator[](size_t index){
                return _data[index];
            }
            const T &operator[](size_t index) const{
                return _data[index];
            }
            size_t size() const{
                return _size;
            }
            void resize(size_t size){
                if(size == _size)
                    return;
                T *new_data = new T[size];
                for(size_t i = 0;i < size;i++){
                    new_data[i] = _data[i];
                }
                delete[] _data;
                _data = new_data;
                _size = size;
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

    //TODO: GuassianBlur
    template<typename View>
    void GuassianBlur(View &view,double sigma){
        
    }
}

#endif // _BTK_PROJECT_PIXUTIL_FLITER_
