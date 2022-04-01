#if !defined(_BTK_PROJECT_PIXUTIL_VIEW_HPP_)
#define _BTK_PROJECT_PIXUTIL_VIEW_HPP_

#include <type_traits>
#include <stdexcept>
#include <cstdint>

#ifndef PIXUTIL_ASSERT
    #define PIXUTIL_ASSERT(x) assert(x)
    #include <cassert>
#endif

#ifdef PIXUTIL_NO_EXCEPTIONS
    #define PIXUTIL_THROW(x) PIXUTIL_ASSERT(false)
#else
    #define PIXUTIL_THROW(x) throw x
#endif

namespace PixUtil{
    //Basic types
    using Uint8 = std::uint8_t;
    using Uint16 = std::uint16_t;
    using Uint32 = std::uint32_t;
    using Uint64 = std::uint64_t;
    using Int8 = std::int8_t;
    using Int16 = std::int16_t;
    using Int32 = std::int32_t;
    using Int64 = std::int64_t;
    //Basic min / max / clamp
    template<typename T>
    constexpr T min(T a,T b) noexcept{
        return a < b ? a : b;
    }
    template<typename T>
    constexpr T max(T a,T b) noexcept{
        return a > b ? a : b;
    }
    template<typename T>
    constexpr T clamp(T value,T min,T max) noexcept{
        return value < min ? min : value > max ? max : value;
    }
    //Color
    #ifndef PIXUTIL_COLOR
    struct Color{
        Uint8 r;
        Uint8 g;
        Uint8 b;
        Uint8 a;
    };
    #else
    using Color = PIXUTIL_COLOR;
    #endif



    template<typename T>
    struct _Pixel{
        T *view;
        int x,y;

        using value_type = typename T::value_type;

        _Pixel() = default;
        _Pixel(T *view,int x,int y):view(view),x(x),y(y){}


        /**
         * @brief assign the value to the pixel
         * 
         * @param pix 
         * @return _Pixel& 
         */
        _Pixel &operator =(value_type pix){
            void *addr = view->index_pixel(
                x,
                y,
                view->bytes_per_pixel()
            );
            view->store_pixel(addr,pix);
            return *this;
        }
        /**
         * @brief Set current pixel by color
         * 
         * @param c 
         * @return _Pixel& 
         */
        _Pixel &operator =(Color c){
            void *addr = view->index_pixel(
                x,
                y,
                view->bytes_per_pixel()
            );
            value_type pix = view->map_color(c);
            view->store_pixel(addr,pix);
            return *this;
        }

        value_type to_pixel() const{
            void *addr = view->index_pixel(
                x,
                y,
                view->bytes_per_pixel()
            );
            return view->load_pixel(addr);
        }
        Color to_color() const{
            value_type pix = to_pixel();
            return view->get_color(pix);
        }

        /**
         * @brief Load the pixel value
         * 
         * @return value_type 
         */
        operator value_type(){
            return to_pixel();
        }
        /**
         * @brief Get Color from pixels
         * 
         * @return Color 
         */
        operator Color(){
            return to_color();
        }
        value_type *operator &() const{
            return  (value_type *)view->index_pixel(
                x,
                y,
                view->bytes_per_pixel()
            );
        }
    };
    template<typename T>
    struct _PixelIterator:private _Pixel<T>{

        _PixelIterator() = default;
        _PixelIterator(T *view,int x,int y):_Pixel<T>(view,x,y){}

        _PixelIterator &operator ++(){
            this->x++;
            return *this;
        }
        _PixelIterator &operator --(){
            this->x--;
            return *this;
        }
        _PixelIterator &operator +=(int n){
            this->x += n;
            return *this;
        }
        _PixelIterator &operator -=(int n){
            this->x -= n;
            return *this;
        }
        _PixelIterator &operator ++(int){
            _PixelIterator tmp = *this;
            this->x++;
            return tmp;
        }
        _PixelIterator &operator --(int){
            _PixelIterator tmp = *this;
            this->x--;
            return tmp;
        }
        _PixelIterator operator +(int n){
            _PixelIterator tmp = *this;
            tmp.x += n;
            return tmp;
        }
        _PixelIterator operator -(int n){
            _PixelIterator tmp = *this;
            tmp.x -= n;
            return tmp;
        }
        bool operator ==(const _PixelIterator &other){
            return this->x == other.x && this->y == other.y;
        }
        bool operator !=(const _PixelIterator &other){
            return !(*this == other);
        }
        //
        _Pixel<T> &operator *(){
            return *this;
        }
        _Pixel<T> *operator ->(){
            return this;
        }
        _Pixel<T> &operator *() const{
            return *this;
        }
        _Pixel<T> *operator ->() const{
            return this;
        }
    };
    /**
     * @brief A view of a row in pixel buffer.
     * 
     * @tparam T 
     */
    template<typename T>
    struct _Row{
        T *view;
        int y;

        _Row() = default;
        _Row(T *view,int y):view(view),y(y){}

        /**
         * @brief Index pixel at x
         * 
         * @param x 
         * @return _Pixel<T> 
         */
        _Pixel<T> operator[](int x) const{
            return _Pixel<T>(view,x,y);
        }
        /**
         * @brief Index a pixel at x(has range check)
         * 
         * @param x 
         * @return _Pixel<T> 
         */
        _Pixel<T> at(int x) const{
            if(!view->has_x(x)){
                throw std::out_of_range("_Row::at() out of range");
            }
            return _Pixel<T>(view,x,y);
        }
        /**
         * @brief Get the begin iterator
         * 
         * @return _PixelIterator<T> 
         */
        _PixelIterator<T> begin() const{
            return _PixelIterator<T>(view,0,y);
        }
        /**
         * @brief Get the end iterator
         * 
         * @return _PixelIterator<T> 
         */
        _PixelIterator<T> end() const{
            return _PixelIterator<T>(view,view->width(),y);
        }
    };
    /**
     * @brief Iterator for a row of pixels
     * 
     * @tparam T 
     */
    template<typename T>
    struct _RowIterator:private _Row<T>{

        _RowIterator() = default;
        _RowIterator(T *view,int y):_Row<T>(view,y){}

        _RowIterator& operator ++(){
            this->y++;
            return *this;
        }
        _RowIterator operator ++(int){
            _RowIterator tmp = *this;
            this->y++;
            return tmp;
        }
        _RowIterator& operator --(){
            this->y--;
            return *this;
        }
        _RowIterator operator --(int){
            _RowIterator tmp = *this;
            this->y--;
            return tmp;
        }

        _RowIterator& operator +=(int n){
            this->y += n;
            return *this;
        }
        _RowIterator& operator -=(int n){
            this->y -= n;
            return *this;
        }
        _RowIterator operator +(int n){
            _RowIterator tmp = *this;
            tmp.y += n;
            return tmp;
        }

        _RowIterator operator -(int n){
            _RowIterator tmp = *this;
            tmp.y -= n;
            return tmp;
        }

        int operator -(const _RowIterator &iter){
            return this->y - iter.y;
        }
        /**
         * @brief Compare two row iterator
         * 
         * @param iter 
         * @return true 
         * @return false 
         */
        bool operator ==(const _RowIterator &iter){
            return this->view == iter.view && this->y == iter.y;
        }
        /**
         * @brief Compare two row iterator
         * 
         * @param iter 
         * @return true 
         * @return false 
         */
        bool operator !=(const _RowIterator &iter){
            return this->view != iter.view || this->y != iter.y;
        }

        _Row<T> &operator *(){
            return *this;
        }
        _Row<T> *operator ->(){
            return this;
        }

        const _Row<T> &operator *() const{
            return *this;
        }
        const _Row<T> *operator ->() const{
            return this;
        }
    };

    /**
     * @brief Generic pixel traits
     * 
     */
    template<typename T>
    struct PixTraits{
        using type = T;

        /**
         * @brief Get the size of the pixel
         * 
         * @return constexpr size_t 
         */
        static constexpr size_t bytes_per_pixel() noexcept{
            return sizeof(T);
        }
        /**
         * @brief Store a pixel at the given address
         * 
         * @param dst 
         * @param src 
         */
        static void store_pixel(void *dst,const T &src){
            *(T*)dst = src;
        }
        /**
         * @brief Load a pixel from the given address
         * 
         * @param src 
         * @return T 
         */
        static T load_pixel(const void *src){
            return *(T*)src;
        }
    };
    /**
     * @brief Traits for RGBA pixels
     * 
     * @tparam T 
     */
    template<typename T>
    struct RGBATraits:public PixTraits<T>{
        static_assert(sizeof(T) == sizeof(Uint32),"RGBATraits::T must be 32 bit");
        //TODO: add endianess support

        static Color get_color(const T pix){
            Color c;
            //Little endian AABBGGRR
            c.r = reinterpret_cast<const Uint8*>(&pix)[0];
            c.g = reinterpret_cast<const Uint8*>(&pix)[1];
            c.b = reinterpret_cast<const Uint8*>(&pix)[2];
            c.a = reinterpret_cast<const Uint8*>(&pix)[3];
            return c;
        }
        static T     map_color(Color c){
            T pix = 0;
            //Little endian AABBGGRR
            reinterpret_cast<Uint8*>(&pix)[0] = c.r;
            reinterpret_cast<Uint8*>(&pix)[1] = c.g;
            reinterpret_cast<Uint8*>(&pix)[2] = c.b;
            reinterpret_cast<Uint8*>(&pix)[3] = c.a;
            return pix;
        }
    };
    /**
     * @brief Traits for RGB pixels
     * 
     * @tparam T 
     */
    template<typename T>
    struct RGBTraits:public PixTraits<T>{
        static_assert(sizeof(T) == sizeof(Uint32),"RGBTraits::T must be 32 bit");
        //TODO: add endianess support

        static Color get_color(const T pix){
            Color c;
            c.r = reinterpret_cast<const Uint8*>(&pix)[0];
            c.g = reinterpret_cast<const Uint8*>(&pix)[1];
            c.b = reinterpret_cast<const Uint8*>(&pix)[2];
            c.a = 255;
            return c;
        }
        static T     map_color(Color c){
            T pix = 0;
            reinterpret_cast<Uint8*>(&pix)[0] = c.r;
            reinterpret_cast<Uint8*>(&pix)[1] = c.g;
            reinterpret_cast<Uint8*>(&pix)[2] = c.b;
            return pix;
        }
    };
    /**
     * @brief View for raw pixels data
     * 
     */
    struct RawView{
        //     ===     BITMAP    ===
        //(0,0)--------Stride----------
        //|
        //|    (x_offset,y_offset)--w--
        //|     |                     |
        //|     h      VIEW AREA      |
        //|     |----------------------

        int x_offset = 0;//< X offset in the bitmap
        int y_offset = 0;//< Y offset in the bitmap
        int stride = 0;//< Stride in the bitmap
        int w = 0;//< Width of the view area
        int h = 0;//< Height of the view area

        void *pixels;
        /**
         * @brief Index a pixel at (x,y)
         * 
         * @param x 
         * @param y 
         * @param bytes_per_pixel 
         * @return void* 
         */
        void *index_pixel(int x,int y,int bytes_per_pixel) const{
            return (Uint8*)pixels + (y + y_offset)* stride + (x + x_offset)* bytes_per_pixel;
        }
        void *index_pixel_at(int x,int y,int bytes_per_pixel) const{
            if (has_point(x,y)){
                return index_pixel(x,y,bytes_per_pixel);
            }
            throw std::out_of_range("RawView::index_pixel_at() out of range");
        }
        bool  has_point(int x,int y) const{
            return x >= 0 && x < w && y >= 0 && y < h;
        }
        bool  has_y(int y) const{
            return y >= 0 && y < h;
        }
        bool  has_x(int x) const{
            return x >= 0 && x < w;
        }
        
    };
    /**
     * @brief Generic pixel view
     * 
     * @tparam Pix 
     * @tparam Traits 
     */
    template<
        typename Pix,
        typename Traits = PixTraits<Pix>
    >
    class View:public RawView,public Traits{
        public:
            View() = default;
            View(const View &) = default;
            View(View &&) = default;
            ~View() = default;
            
            View &operator=(const View &) = default;
            View &operator=(View &&) = default;
        public:
            /**
             * @brief Construct a new View object
             * 
             * @param pixels 
             * @param x_offset 
             * @param y_offset 
             * @param stride 
             * @param w 
             * @param h 
             * @param traits 
             */
            View(void *pixels,int x_offset,int y_offset,int stride,int w,int h,Traits traits = Traits())
                :Traits(traits){
                
                this->pixels = pixels;
                this->x_offset = x_offset;
                this->y_offset = y_offset;
                this->stride = stride;
                this->w = w;
                this->h = h;
            }
            /**
             * @brief Construct a new View object with w h and stride
             * 
             * @param pixels 
             * @param w 
             * @param h 
             * @param stride 
             * @param traits 
             */
            View(void *pixels,int w,int h,int stride,Traits traits = Traits())
                :Traits(traits){
                
                this->pixels = pixels;
                this->x_offset = 0;
                this->y_offset = 0;
                this->stride = stride;
                this->w = w;
                this->h = h;
            }
            /**
             * @brief Construct a new View object
             * 
             * @param pixels 
             * @param w 
             * @param h 
             * @param traits 
             */
            View(void *pixels,int w,int h,Traits traits = Traits())
                :Traits(traits){
                
                this->pixels = pixels;
                this->w = w;
                this->h = h;
                this->stride = w * Traits::bytes_per_pixel();
            }
        public:
            View subview(int x,int y,int w,int h) const{
                x = clamp(x,0,this->w);
                y = clamp(y,0,this->h);
                w = clamp(w,0,this->w - x);
                h = clamp(h,0,this->h - y);

                return View(
                    pixels,
                    x_offset + x,
                    y_offset + y,
                    stride,
                    w,
                    h,
                    *this
                );
            }
        public:
            int width() const{
                return w;
            }
            int height() const{
                return h;
            }
            const Traits &traits() const{
                return *this;
            }
        public:
            using value_type = Pix;
            using iterator = _RowIterator<View>;
            using const_iterator = _RowIterator<const View>;

            using reference = _Row<View>;
            using const_reference = _Row<const View>;

            using pixel_reference = _Pixel<View>;
            using const_pixel_reference = _Pixel<const View>;

            iterator begin(){
                return iterator(this,0);
            }
            iterator end(){
                return iterator(this,h);
            }
            const_iterator begin() const{
                return const_iterator(this,0);
            }
            const_iterator end() const{
                return const_iterator(this,h);
            }

            reference operator[](int y){
                return reference(this,y);
            }
            const_reference operator[](int y) const{
                return const_reference(this,y);
            }
            /**
             * @brief Index a row by giving y
             * 
             * @param y 
             * @return reference 
             */
            reference at(int y){
                //Check has y
                if(!has_y(y)){
                    throw std::out_of_range("View::at() out of range");
                }
                return reference(this,y);
            }
            const_reference at(int y) const{
                //Check has y
                if(!has_y(y)){
                    throw std::out_of_range("View::at() out of range");
                }
                return const_reference(this,y);
            }
            /**
             * @brief Index a pixel by giving position
             * 
             * @param x 
             * @param y 
             * @return pixel_reference 
             */
            pixel_reference at(int x,int y){
                if(!has_point(x,y)){
                    throw std::out_of_range("View::at() out of range");
                }
                return pixel_reference(this,x,y);
            }
            const_pixel_reference at(int x,int y) const{
                if(!has_point(x,y)){
                    throw std::out_of_range("View::at() out of range");
                }
                return const_pixel_reference(this,x,y);
            }
            /**
             * @brief Index a pixel by giving position(unchecked)
             * 
             * @param x 
             * @param y 
             * @return pixel_reference 
             */
            pixel_reference unsafe_at(int x,int y){
                return pixel_reference(this,x,y);
            }
            const_pixel_reference unsafe_at(int x,int y) const{
                return const_pixel_reference(this,x,y);
            }

    };

    //RGBA View
    using RGBAView = View<Uint32,RGBATraits<Uint32>>;
    using RGBView = View<Uint32,RGBTraits<Uint32>>;

    /**
     * @brief Copy pixels from a view to another view(must have same size)
     * 
     * @tparam View1 
     * @tparam View2 
     * @param dst 
     * @param src 
     */
    template<typename View1,typename View2>
    void CopyPixels(View1 &dst,const View2 &src){
        PIXUTIL_ASSERT(dst.width() == src.width());
        PIXUTIL_ASSERT(dst.height() == src.height());

        for(int y = 0;y < src.height();++ y){
            for(int x = 0;x < src.width(); ++x){
                dst[y][x] = src[y][x].to_pixel();
            }
        }
    }
    /**
     * @brief Fill a rectangle with a color / pixel
     * 
     * @tparam View 
     * @tparam T 
     * @param view 
     * @param x 
     * @param y 
     * @param w 
     * @param h 
     * @param p 
     */
    template<typename View,typename T>
    void FillRect(View &view,int x,int y,int w,int h,T p){
        x = clamp(x,0,view.width());
        y = clamp(y,0,view.height());
        w = clamp(w,0,view.width() - x);
        h = clamp(h,0,view.height() - y);

        for(int yy = y;yy < y + h;++ yy){
            for(int xx = x;xx < x + w;++ xx){
                view[yy][xx] = p;
            }
        }
    }
    /**
     * @brief Namespace for blend (dst,src) => new
     * 
     */
    namespace Blend{
        /**
         * @brief None blend (dst,src) => src
         * 
         * @param src 
         * @return src 
         */
        inline Color None(Color,Color src) noexcept{
            return src;
        }
    } // namespace Blend
    
} // namespace PixUtil


#if defined(PIXUTIL_SDL_EXTERNAL)

namespace PixUtil{
    /**
     * @brief Traits for manage surface
     * 
     */
    class SDLSurfaceTraits{
        public:
            SDLSurfaceTraits() = default;
            SDLSurfaceTraits(SDL_Surface *surf){
                surface = surf;
                if(surface != nullptr){
                    surface->refcount += 1;
                }
            }
            SDLSurfaceTraits(const SDLSurfaceTraits &t){
                surface = t.surface;
                if(surface != nullptr){
                    surface->refcount += 1;
                }
            }
            SDLSurfaceTraits(SDLSurfaceTraits &&t){
                surface = t.surface;
                t.surface = nullptr;
            }
            ~SDLSurfaceTraits(){
                SDL_FreeSurface(surface);
            }

            //Operator
            SDLSurfaceTraits &operator=(SDLSurfaceTraits &&t){
                if(this != &t){
                    SDL_FreeSurface(surface);
                    surface = t.surface;
                    t.surface = nullptr;
                }
                return *this;
            }
            SDLSurfaceTraits &operator=(const SDLSurfaceTraits &t){
                if(this != &t){
                    SDL_FreeSurface(surface);
                    surface = t.surface;
                    if(surface != nullptr){
                        surface->refcount += 1;
                    }
                }
                return *this;
            }
        public:
            size_t bytes_per_pixel() const{
                return surface->format->BytesPerPixel;
            }
            /**
             * @brief Load pixels by BytesPerPixel
             * 
             * @param pixel 
             * @return Uint32 
             */
            Uint32 load_pixel(const void *pixel) const{
                Uint32 pix;
                switch(bytes_per_pixel()){
                    case 1:
                        pix = *(Uint8 *)pixel;
                        break;
                    case 2:
                        pix = *(Uint16 *)pixel;
                        break;
                    case 3:
                        pix = *(Uint32 *)pixel;
                        break;
                    case 4:
                        pix = *(Uint32 *)pixel;
                        break;
                    default:
                        pix = 0;
                        break;
                }
                return pix;
            }
            void   store_pixel(void *pixel,Uint32 p) const{
                switch(bytes_per_pixel()){
                    case 1:
                        *(Uint8 *)pixel = p;
                        break;
                    case 2:
                        *(Uint16 *)pixel = p;
                        break;
                    case 3:
                        *(Uint32 *)pixel = p;
                        break;
                    case 4:
                        *(Uint32 *)pixel = p;
                        break;
                    default:
                        break;
                }
            }
            Color get_color(Uint32 pix){
                Color c;
                SDL_GetRGBA(
                    pix,
                    surface->format,
                    &c.r,
                    &c.g,
                    &c.b,
                    &c.a
                );
                return c;
            }
            Uint32 map_color(Color c){
                return SDL_MapRGBA(
                    surface->format,
                    c.r,
                    c.g,
                    c.b,
                    c.a
                );
            }
        private:
            SDL_Surface *surface;
    };
    /**
     * @brief Wrapper for SDL_Surface and View
     * 
     */
    class SDLSurfaceView:public View<Uint32,SDLSurfaceTraits>{
        public:
            using ViewBase = View<Uint32,SDLSurfaceTraits>;
            using ViewBase::View;
            /**
             * @brief Construct a new SDLSurfaceView object
             * 
             * @param surf The surface pointer(nullptr is not allowed)
             */
            SDLSurfaceView(SDL_Surface *surf):
                View(surf->pixels,surf->w,surf->h,surf->pitch,SDLSurfaceTraits(surf)){

            }
            SDLSurfaceView(const ViewBase &b):ViewBase(b){}
            SDLSurfaceView(ViewBase &&b):ViewBase(b){}
            ~SDLSurfaceView() = default;

            SDLSurfaceView subview(int x,int y,int w,int h) const{
                return ViewBase::subview(x,y,w,h);
            }
    };
} // namespace PixUtil


#endif // PIXUTIL_SDL_EXTERNAL

#endif // _BTK_PROJECT_PIXUTIL_VIEW_HPP_
