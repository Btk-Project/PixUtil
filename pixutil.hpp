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
    /**
     * @brief C++ 11 void_t
     * 
     */
    template<typename...>
    using void_t = void;
    /**
     * @brief Pixel reference
     * 
     * @tparam T 
     */
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
     * @brief Internal ShiftTraits
     * 
     * @tparam Mask 
     * @tparam Shift 
     * @tparam !(Mask & 0x01) 
     */
    template<Uint32 Mask,Uint8 Shift = 0,bool v = !(Mask & 0x01) && (Mask != 0)>
    struct _ShiftTraits;

    template<Uint32 Mask,Uint8 Shift>
    struct _ShiftTraits<Mask,Shift,false>{
        static constexpr Uint32  mask_value = Mask;
        static constexpr Uint8  shift_value = Shift;
    };
    template<Uint32 Mask,Uint8 Shift>
    struct _ShiftTraits<Mask,Shift,true>{
        static constexpr Uint32  mask_value = _ShiftTraits<(Mask >> 1),Shift + 1>::mask_value;
        static constexpr Uint8  shift_value = _ShiftTraits<(Mask >> 1),Shift + 1>::shift_value;
    };
    /**
     * @brief Internal Loss Traits
     * 
     * @tparam Mask 
     * @tparam Loss 
     * @tparam (Mask & 0x01) 
     */
    template<Uint32 Mask,Uint8 Loss = 8,bool v = (Mask & 0x01) && (Mask != 0)>
    struct _LossTraits;

    template<Uint32 Mask,Uint8 Loss>
    struct _LossTraits<Mask,Loss,false>{
        static constexpr Uint32 mask_value = Mask;
        static constexpr Uint8  loss_value = Loss;
    };
    template<Uint32 Mask,Uint8 Loss>
    struct _LossTraits<Mask,Loss,true>{
        static constexpr Uint32 mask_value = _LossTraits<(Mask >> 1),Loss - 1>::mask_value;
        static constexpr Uint8  loss_value = _LossTraits<(Mask >> 1),Loss - 1>::loss_value;
    };
    /**
     * @brief Traits for Get Shift from a Uint32 Mask
     * 
     * @tparam Mask 
     */
    template<Uint32 Mask>
    using ShiftTraits = _ShiftTraits<Mask>;
    /**
     * @brief Traits for Get Loss from a Uint32 Mask
     * 
     * @tparam Mask 
     */
    template<Uint32 Mask>
    using LossTraits = _LossTraits<_ShiftTraits<Mask>::mask_value>;
    /**
     * @brief Traits for Generic RGBA Pixel like SDL_PixelFormat
     * 
     * @tparam T Pixel type
     * @tparam BitsPerPixel
     * @tparam Rmask
     * @tparam Gmask
     * @tparam Bmask
     * @tparam Amask
     */
    template<
        typename T,//< Pixel type
        size_t _BitsPerPixel,
        Uint32 _Rmask,
        Uint32 _Gmask,
        Uint32 _Bmask,
        Uint32 _Amask
    >
    struct RGBACommonTraits{

        //Get bytes per pixel
        static constexpr size_t BytesPerPixel = (_BitsPerPixel + 7) / 8;
        static constexpr size_t BitsPerPixel = _BitsPerPixel;
        //Get Rmask,Gmask,Bmask,Amask
        static constexpr Uint32 Rmask = _Rmask;
        static constexpr Uint32 Gmask = _Gmask;
        static constexpr Uint32 Bmask = _Bmask;
        static constexpr Uint32 Amask = _Amask;
        //Get Rloss,Gloss,Bloss,Aloss
        static constexpr Uint8  Rloss = LossTraits<Rmask>::loss_value;
        static constexpr Uint8  Gloss = LossTraits<Gmask>::loss_value;
        static constexpr Uint8  Bloss = LossTraits<Bmask>::loss_value;
        static constexpr Uint8  Aloss = LossTraits<Amask>::loss_value;
        //Get Rshift,Gshift,Bshift,Ashift
        static constexpr Uint8  Rshift = ShiftTraits<Rmask>::shift_value;
        static constexpr Uint8  Gshift = ShiftTraits<Gmask>::shift_value;
        static constexpr Uint8  Bshift = ShiftTraits<Bmask>::shift_value;
        static constexpr Uint8  Ashift = ShiftTraits<Amask>::shift_value;

        static_assert(sizeof(T) >= BytesPerPixel,"Invalid BytesPerPixel");
        static_assert(sizeof(T) == sizeof(Uint32),"Invalid BytesPerPixel");

        //Load / save pixel
        static constexpr size_t bytes_per_pixel() noexcept{
            return BytesPerPixel;
        }
        static void store_pixel(void *dst,const T &src){
            switch(BytesPerPixel){
                case 1:
                    *(Uint8*)dst = src;
                    break;
                case 2:
                    *(Uint16*)dst = src;
                    break;
                case 3:
                    *(Uint8*)dst = src;
                    *(Uint8*)((Uint8*)dst + 1) = src >> 8;
                    *(Uint8*)((Uint8*)dst + 2) = src >> 16;
                    break;
                case 4:
                    *(Uint32*)dst = src;
                    break;
                default:
                    break;
            }
        }
        static T    load_pixel(const void *src){
            T pix;
            switch(BytesPerPixel){
                case 1:
                    pix = *(T*)src;
                    break;
                case 2:
                    pix = *(T*)src;
                    break;
                case 3:
                    pix = *(T*)src;
                    break;
                case 4:
                    pix = *(T*)src;
                    break;
                default:
                    break;
            }
            return pix;
        }

        static Color get_color(const T pix){
            Color c;
            c.r = (pix & Rmask) >> Rshift;
            c.g = (pix & Gmask) >> Gshift;
            c.b = (pix & Bmask) >> Bshift;
            c.a = (pix & Amask) >> Ashift;

            //Process with Loss
            c.r = (c.r << Rloss) | (c.r >> (8 - Rloss));
            c.g = (c.g << Gloss) | (c.g >> (8 - Gloss));
            c.b = (c.b << Bloss) | (c.b >> (8 - Bloss));
            c.a = (c.a << Aloss) | (c.a >> (8 - Aloss));
            return c;
        }
        static T map_color(const Color &c){
            T pix;
            //Map Color by Loss Shift and Mask
            pix =  (c.r >> Rloss) << Rshift;
            pix |= (c.g >> Gloss) << Gshift;
            pix |= (c.b >> Bloss) << Bshift;
            //Because Alpha may not be used,so we need to process it by Amask
            pix |= Uint32(c.a >> Aloss) << Ashift & Amask;

            return pix;
        }
    };
    /**
     * @brief Traits for RGBA32 pixels
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
     * @brief Traits for RGB24 pixels (24 bit)
     * 
     * @tparam T 
     */
    template<typename T>
    struct RGBTraits{
        static_assert(sizeof(T) == sizeof(Uint32),"RGBTraits::T must be 32 bit");
        //TODO: add endianess support

        static constexpr size_t bytes_per_pixel() noexcept{
            return 3;
        }
        static void store_pixel(void *dst,const T &src){
            reinterpret_cast<Uint8*>(dst)[0] = reinterpret_cast<const Uint8*>(&src)[0];
            reinterpret_cast<Uint8*>(dst)[1] = reinterpret_cast<const Uint8*>(&src)[1];
            reinterpret_cast<Uint8*>(dst)[2] = reinterpret_cast<const Uint8*>(&src)[2];
        }
        static T load_pixel(const void *src){
            T pix = 0;
            reinterpret_cast<Uint8*>(&pix)[0] = reinterpret_cast<const Uint8*>(src)[0];
            reinterpret_cast<Uint8*>(&pix)[1] = reinterpret_cast<const Uint8*>(src)[1];
            reinterpret_cast<Uint8*>(&pix)[2] = reinterpret_cast<const Uint8*>(src)[2];
            return pix;
        }

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
     * @brief Traits for Grayscale pixels (8 bit)
     * @note The get / map is only for erase the pixels type(donot use it for conversion)
     * @tparam T 
     */
    template<typename T>
    struct GrayTraits:public PixTraits<T>{
        static_assert(sizeof(T) == sizeof(Uint8),"GrayTraits::T must be 8 bit");

        static Color get_color(const T pix){
            //Convert grayscale to RGB
            Color c;
            c.r = c.g = c.b = pix;
            c.a = 255;
            return c;
        }
        static T     map_color(Color c){
            //Convert RGBA to grayscale
            return (Uint32(c.r) + Uint32(c.g) + Uint32(c.b)) / 3;
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
    //RGBA Common View
    template<size_t Bits,Uint32 Rmask,Uint32 Gmask,Uint32 Bmask,Uint32 Amask = 0>
    using RGBACommonView = View<Uint32,RGBACommonTraits<Uint32,Bits,Rmask,Gmask,Bmask,Amask>>;

    //RGBA View
    using RGBAView = View<Uint32,RGBATraits<Uint32>>;
    using GrayView = View<Uint8,GrayTraits<Uint8>>;
    using RGBView = View<Uint32,RGBTraits<Uint32>>;

    #ifdef PIXUTIL_COMMON_PIXELFORMAT
    //Add common pixel format by using Rmask,Gmask,Bmask,Amask
    using BGRAView = RGBACommonView<32,0xff0000,0x00ff00,0x0000ff,0xff000000>;
    // using RGBAView = RGBACommonView<32,0x000000ff,0x0000ff00,0x00ff0000,0xff000000>;
    using ARGBView = RGBACommonView<32,0x00ff0000,0x0000ff00,0x000000ff,0xff000000>;
    using ABGRView = RGBACommonView<32,0x0000ff00,0x00ff0000,0xff000000,0x000000ff>;
    using XBGRView = RGBACommonView<32,0x000000ff,0x0000ff00,0x00ff0000,0x00000000>;
    using XRGBView = RGBACommonView<32,0x00ff0000,0x0000ff00,0x000000ff,0x00000000>;
    using BGRXView = RGBACommonView<32,0x00ff0000,0x0000ff00,0x000000ff,0x00000000>;
    using RGBXView = RGBACommonView<32,0x000000ff,0x0000ff00,0x00ff0000,0x00000000>;
    //24
    using BGRView = RGBACommonView<24,0x00ff0000,0x0000ff00,0x000000ff,0x00000000>;
    // using RGBView = RGBACommonView<24,0x000000ff,0x0000ff00,0x00ff0000,0x00000000>;
    //16
    using RGB565View = RGBACommonView<16,0x0000f800,0x000007e0,0x0000001f,0x00000000>;
    using BGR565View = RGBACommonView<16,0x0000001f,0x000007e0,0x0000f800,0x00000000>;
    using RGB555View = RGBACommonView<16,0x00007c00,0x000003e0,0x0000001f,0x00000000>;
    using BGR555View = RGBACommonView<16,0x0000001f,0x000003e0,0x00007c00,0x00000000>;
    //8
    using RGB332View = RGBACommonView<8,0x000000e0,0x0000001c,0x00000003,0x00000000>;
    //4
    using RGB444View = RGBACommonView<4,0x00000f00,0x000000f0,0x0000000f,0x00000000>;
    using BGR444View = RGBACommonView<4,0x0000000f,0x000000f0,0x00000f00,0x00000000>;
    using RGB422View = RGBACommonView<4,0x00000f00,0x000000f0,0x0000000f,0x00000000>;
    using BGR422View = RGBACommonView<4,0x0000000f,0x000000f0,0x00000f00,0x00000000>;
    //2
    using RGB211View = RGBACommonView<2,0x00000c00,0x00000300,0x00000001,0x00000000>;
    using BGR211View = RGBACommonView<2,0x00000001,0x00000300,0x00000c00,0x00000000>;
    #endif

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
     * @brief Convert RGB to Gray
     * 
     * @tparam View1 
     * @tparam View2 
     * @param dst 
     * @param src 
     */
    template<typename View1,typename View2>
    void RGBToGray(View1 &dst,const View2 &src){
        PIXUTIL_ASSERT(dst.width() == src.width());
        PIXUTIL_ASSERT(dst.height() == src.height());

        for(int y = 0;y < src.height();++ y){
            for(int x = 0;x < src.width(); ++x){
                Color c = src[y][x];
                Uint8 pix = c.r * 0.299 + c.g * 0.587 + c.b * 0.114;
                dst[y][x] = pix;
            }
        }
    }
    /**
     * @brief Using Nearest sampling to scale a view to another view
     * 
     * @tparam View1 
     * @tparam View2 
     * @param dst 
     * @param src 
     */
    template<typename View1,typename View2>
    void NearestScale(View1 &dst,const View2 &src){
        //Nearest sampling
        int out_w = dst.width();
        int out_h = dst.height();
        int in_w = src.width();
        int in_h = src.height();
        //Begin sampling
        for(int y = 0;y < out_h;++ y){
            for(int x = 0;x < out_w;++ x){
                //Get the source pixel
                int in_x = x * in_w / out_w + 0.5;
                int in_y = y * in_h / out_h + 0.5;

                //Assert the pixel is in the source view
                PIXUTIL_ASSERT(in_x >= 0 && in_x < in_w);
                PIXUTIL_ASSERT(in_y >= 0 && in_y < in_h);
                
                dst[y][x] = src[in_y][in_x].to_color();
            }
        }
    }
    /**
     * @brief Using Bilinear sampling to scale a view to another view
     * 
     * @tparam View1 
     * @tparam View2 
     * @param dst 
     * @param src 
     */
    template<typename View1,typename View2>
    void BilinearScale(View1 &dst,const View2 &src){
        //Bilinear sampling
        int out_w = dst.width();
        int out_h = dst.height();
        int in_w = src.width();
        int in_h = src.height();
        //Begin sampling
        for(int y = 0;y < out_h;++ y){
            for(int x = 0;x < out_w;++ x){
                //Get the source pixel
                int in_x = x * in_w / out_w;
                int in_y = y * in_h / out_h;

                //Assert the pixel is in the source view
                PIXUTIL_ASSERT(in_x >= 0 && in_x < in_w);
                PIXUTIL_ASSERT(in_y >= 0 && in_y < in_h);

                //Get the source pixel
                Color c00 = src[in_y][in_x];
                Color c01 = src[in_y][in_x + 1];
                Color c10 = src[in_y + 1][in_x];
                Color c11 = src[in_y + 1][in_x + 1];

                //Get the weight
                float wx = x * in_w / out_w - in_x;
                float wy = y * in_h / out_h - in_y;

                //Get the color
                Color c;
                c.r = c00.r * (1 - wx) * (1 - wy) + c01.r * wx * (1 - wy) + c10.r * (1 - wx) * wy + c11.r * wx * wy;
                c.g = c00.g * (1 - wx) * (1 - wy) + c01.g * wx * (1 - wy) + c10.g * (1 - wx) * wy + c11.g * wx * wy;
                c.b = c00.b * (1 - wx) * (1 - wy) + c01.b * wx * (1 - wy) + c10.b * (1 - wx) * wy + c11.b * wx * wy;
                c.a = c00.a * (1 - wx) * (1 - wy) + c01.a * wx * (1 - wy) + c10.a * (1 - wx) * wy + c11.a * wx * wy;
                
                dst[y][x] = c;
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

#include <SDL2/SDL_surface.h>

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
            //FIXME :valgrind reported using uninitialised value here
            Color get_color(Uint32 pix) const{
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
            Uint32 map_color(Color c) const{
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

    //Helper for show view in a window,like opencv

    template<class View>
    void SDLShowView(const View &view,const char *title = nullptr){
        static int n = 0;
        if(title == nullptr){
            char tmp_buffer[128] = {0};
            size_t len;
            SDL_snprintf(
                tmp_buffer,
                sizeof(tmp_buffer),
                "SDLView [%d] - %d * %d",
                n,
                view.width(),
                view.height()
            );
            n += 1;
            //Copy to title
            len = SDL_strlen(tmp_buffer);
            title = SDL_stack_alloc(char,len + 1);
            SDL_strlcpy(const_cast<char*>(title),tmp_buffer,len + 1);
        }
        //Create a surface and copy the view to it
        SDL_Surface *surface = SDL_CreateRGBSurfaceWithFormat(
            0,
            view.width(),
            view.height(),
            32,
            SDL_PIXELFORMAT_RGBA32
        );
        SDLSurfaceView sdl_view(surface);
        //Convert
        for(int y = 0;y < view.height();++y){
            for(int x = 0;x < view.width();++x){
                sdl_view[y][x] = view[y][x].to_color();
            }
        }
        //Create done
        SDL_Window *window = SDL_CreateWindow(
            title,
            SDL_WINDOWPOS_UNDEFINED,
            SDL_WINDOWPOS_UNDEFINED,
            view.width(),
            view.height(),
            SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE
        );
        SDL_Renderer *renderer = SDL_CreateRenderer(
            window,
            -1,
            0
        );
        SDL_Texture *texture = SDL_CreateTextureFromSurface(
            renderer,
            surface
        );
        //Bind
        SDL_SetWindowData(window,"SDLViewTexture",texture);
        //End
        SDL_FreeSurface(surface);
    }
    /**
     * @brief Wait for all window closed
     * 
     */
    inline void SDLWaitForQuit(){
        SDL_Event event;
        while(SDL_WaitEvent(&event)){
            if(event.type == SDL_QUIT){
                break;
            }
            else if(event.type == SDL_WINDOWEVENT){
                SDL_Window *window = SDL_GetWindowFromID(event.window.windowID);
                if(window == nullptr){
                    continue;
                }
                if(event.window.event == SDL_WINDOWEVENT_EXPOSED){
                    //Get texture
                    SDL_Renderer *renderer = SDL_GetRenderer(window);
                    SDL_Texture  *texture = static_cast<SDL_Texture *>(
                        SDL_GetWindowData(window,"SDLViewTexture")
                    );
                    int win_w,win_h;
                    int tex_w,tex_h;

                    SDL_QueryTexture(texture,nullptr,nullptr,&tex_w,&tex_h);
                    SDL_GetWindowSize(window,&win_w,&win_h);
                    //Calc the texture rect to fit the window
                    //Keep the aspect ratio
                    SDL_Rect rect;
                    if(tex_w * win_h > tex_h * win_w){
                        rect.w = win_w;
                        rect.h = tex_h * win_w / tex_w;
                    }else{
                        rect.w = tex_w * win_h / tex_h;
                        rect.h = win_h;
                    }
                    rect.x = (win_w - rect.w) / 2;
                    rect.y = (win_h - rect.h) / 2;
                    //Draw 
                    SDL_RenderClear(renderer);
                    SDL_RenderCopy(renderer,texture,nullptr,&rect);
                    SDL_RenderPresent(renderer);
                    continue;
                }
                else if(event.window.event == SDL_WINDOWEVENT_CLOSE){
                    //Cleanup
                    SDL_DestroyTexture(static_cast<SDL_Texture *>(
                        SDL_GetWindowData(window,"SDLViewTexture")
                    ));
                    SDL_DestroyRenderer(
                        SDL_GetRenderer(window)
                    );
                    SDL_DestroyWindow(window);
                    continue;
                }
            }
        }
    }
} // namespace PixUtil


#endif // PIXUTIL_SDL_EXTERNAL

#endif // _BTK_PROJECT_PIXUTIL_VIEW_HPP_
