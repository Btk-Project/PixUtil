#if !defined(_BTK_PROJECT_PIXUTIL_IMAGE_HPP_)
#define _BTK_PROJECT_PIXUTIL_IMAGE_HPP_
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <cstddef>
#include <cstdio>

//All image load entry macro
#define _PIXUTIL_LIBPNG_LOAD_ENTRY
#define _PIXUTIL_WIC_LOAD_ENTRY
#define _PIXUTIL_STB_LOAD_ENTRY
//Does has_include supported?
#ifdef __has_include
    #define PIXUTIL_HAS_INCLUDE(X) __has_include(X)
#else
    #define PIXUTIL_HAS_INCLUDE(X) false
#endif
//Assert
#ifndef PIXUTIL_ASSERT
    #define PIXUTIL_ASSERT(X) assert(X)
    #include <cassert>
#endif
//Malloc / Free
#ifndef PIXUTIL_IMG_MALLOC
    #define PIXUTIL_IMG_REALLOC std::realloc
    #define PIXUTIL_IMG_MALLOC std::malloc
    #define PIXUTIL_IMG_FREE std::free
#endif
//Memset
#ifndef PIXUTIL_IMG_MEMSET
    #define PIXUTIL_IMG_MEMSET ::memset
#endif
//Wincodec
#ifdef _WIN32
    #ifdef _MSC_VER
        #pragma comment(lib,"ole32.lib")
    #endif
    #define NOMINMAX
    #include <wrl/client.h>
    #include <wincodec.h>
#endif

namespace PixImage{
    inline void *Malloc(size_t size){
        return PIXUTIL_IMG_MALLOC(size);
    }
    inline void *Realloc(void *ptr, size_t size){
        return PIXUTIL_IMG_REALLOC(ptr, size);
    }
    inline void Free(void *ptr){
        PIXUTIL_IMG_FREE(ptr);
    }
    using ssize_t = std::ptrdiff_t;
    /**
     * @brief Simplest image format
     * 
     */
    namespace Format{
        enum {
            Gray  = 1,
            GrayA = 2,
            RGB   = 3,
            RGBA  = 4
        };
    }
    struct Bitmap{
        int w = 0;
        int h = 0;
        int pitch = 0;
        int channels = 0;
        void *pixels = nullptr;

        bool empty() const{
            return w == 0 || h == 0 || pixels == nullptr;
        }
    };
    struct Animation{
        int     w = 0;
        int     h = 0;
        int     nframes = 0;
        int    *durations = nullptr;
        Bitmap *bitmaps = nullptr;
    };
    struct IOCallback{
        void *uptr;
        ssize_t (*read)(void *uptr, void *ptr, size_t size);
        ssize_t (*write)(void *uptr, const void *ptr, size_t size);
        ssize_t (*seek)(void *uptr, ssize_t offset, int whence);
        ssize_t (*tell)(void *uptr);
        void    (*close)(void *uptr);
    };
    /**
     * @brief Alloc a new animation
     * 
     * @param w 
     * @param h 
     * @param frames 
     * @return Animation* 
     */
    inline Animation *AllocAnimation(int w,int h,int channals,int frames){
        Animation *anim = (Animation*)Malloc(sizeof(Animation));
        anim->w = w;
        anim->h = h;
        anim->nframes = frames;
        anim->durations = (int*)Malloc(sizeof(int) * frames);
        anim->bitmaps = (Bitmap*)Malloc(sizeof(Bitmap) * frames);
        //Zero memory
        for(int i = 0;i < frames;++i){
            anim->durations[i] = 0;
            anim->bitmaps[i].w = w;
            anim->bitmaps[i].h = h;
            anim->bitmaps[i].pitch = w * channals;
            anim->bitmaps[i].channels = channals;
            anim->bitmaps[i].pixels = Malloc(w * h * channals);
        }
        return anim;
    }
    inline Animation *ReallocAnimation(Animation *anim,int frames){
        if(anim == nullptr){
            return nullptr;
        }
        int old_frames = anim->nframes;
        anim->nframes = frames;
        anim->durations = (int*)Realloc(anim->durations, sizeof(int) * frames);
        anim->bitmaps = (Bitmap*)Realloc(anim->bitmaps, sizeof(Bitmap) * frames);
        //TODO : prepare for new frames
        return anim;
    }
    /**
     * @brief Free a animation
     * 
     * @param anim 
     * @param free_bitmaps did you want to free bitmaps?
     */
    inline void      FreeAnimation(Animation *anim,bool free_bitmaps = true){
        if(free_bitmaps){
            for(int i = 0;i < anim->nframes;++i){
                Free(anim->bitmaps[i].pixels);
            }
        }
        Free(anim->durations);
        Free(anim->bitmaps);
        Free(anim);
    }
    
    inline void     *AllocFrame(int w,int h,int channels){
        return Malloc(w * h * channels);
    }

    inline void      FreeBitmap(Bitmap bitmap){
        Free(bitmap.pixels);
    }

    //Main functions
    inline Bitmap     LoadFromCallback(const IOCallback &cb,int want = Format::RGBA);
    inline Bitmap     LoadFromFilename(const char *filename,int want = Format::RGBA);
    inline Bitmap     LoadFromMemory(const void *data,size_t size,int want = Format::RGBA);
    inline Animation* LoadAnimation(const IOCallback &cb,int want = Format::RGBA);
    inline IOCallback WrapFILE(FILE *fp);
    inline IOCallback WrapMem(const void *mem,size_t n);
    inline IOCallback WrapMem(void *mem,size_t n);
} // namespace PixImage

#ifdef _WIN32
namespace PixImage{
namespace Wincodec{
    #ifndef PIXUTIL_WIC_NO_BUILTIN_GUID
    constexpr GUID ID_Factory = {0x317d06e8, 0x5f24, 0x433d, 0xbd, 0xf7, 0x79, 0xce,0x68, 0xd8, 0xab, 0xc2};
    constexpr GUID ID_Gray  = {0x6fddc324, 0x4e03, 0x4bfe, 0xb1, 0x85, 0x3d, 0x77, 0x76, 0x8d, 0xc9, 0x08};
    constexpr GUID ID_GrayA = {0x6fddc324, 0x4e03, 0x4bfe, 0xb1, 0x85, 0x3d, 0x77, 0x76, 0x8d, 0xc9, 0x0b};
    constexpr GUID ID_RGB   = {0x6fddc324, 0x4e03, 0x4bfe, 0xb1, 0x85, 0x3d, 0x77, 0x76, 0x8d, 0xc9, 0x0d};
    constexpr GUID ID_RGBA  = {0x3cc4a650, 0xa527, 0x4d37, 0xa9, 0x16, 0x31, 0x42, 0xc7, 0xeb, 0xed, 0xba};
    #else
    constexpr auto ID_Factory = CLSID_WICImagingFactory;
    constexpr auto ID_Gray  = GUID_WICPixelFormat8bppGray;
    constexpr auto ID_GrayA = GUID_WICPixelFormat16bppGray;
    constexpr auto ID_RGB   = GUID_WICPixelFormat24bppRGB;
    constexpr auto ID_RGBA  = GUID_WICPixelFormat32bppPRGBA;
    #endif
    template<typename T>
    using ComPtr = Microsoft::WRL::ComPtr<T>;
    /**
     * @brief Get the Factory object
     * 
     * @return IWICImagingFactory* 
     */
    inline IWICImagingFactory *GetFactory(){
        struct Holder{
            Holder(){
                HRESULT hr;
                hr = CoCreateInstance(
                    ID_Factory,
                    nullptr,
                    CLSCTX_INPROC_SERVER,
                    IID_PPV_ARGS(&factory)
                );
                //Check is com initialized
                if(FAILED(hr)){
                    //Try initialize com
                    hr = CoInitializeEx(nullptr, COINIT_MULTITHREADED);
                    if(FAILED(hr)){
                        return;
                    }
                    //Create factory
                    hr = CoCreateInstance(
                        ID_Factory,
                        nullptr,
                        CLSCTX_INPROC_SERVER,
                        IID_PPV_ARGS(&factory)
                    );
                    if(FAILED(hr)){
                        return;
                    }
                    //< Ask destructor to release com
                    co_initialized = true;
                }
            }
            ~Holder(){
                //Check has Wincodec module

                if(factory != nullptr){
                    if(::GetModuleHandleA("WindowsCodecs.dll") == nullptr){
                        return;
                    }
                    factory->Release();
                }
                if(co_initialized){
                    CoUninitialize();
                }
            }

            IWICImagingFactory *factory = nullptr;
            bool co_initialized = false;
        };
        static thread_local Holder holder;
        return holder.factory;
    }
    inline IStream *CallbackToIStream(const IOCallback &cb){
        cb.seek(cb.uptr,0,SEEK_END);
        size_t n = cb.tell(cb.uptr);
        cb.seek(cb.uptr,0,SEEK_SET);

        //Alloc memory 
        HGLOBAL hmem = GlobalAlloc(GMEM_MOVEABLE,n);
        if(hmem == nullptr){
            return nullptr;
        }
        void *ptr = GlobalLock(hmem);
        if(ptr == nullptr){
            GlobalFree(hmem);
            return nullptr;
        }
        //Read data
        if(cb.read(cb.uptr,ptr,n) != n){
            GlobalUnlock(hmem);
            GlobalFree(hmem);
            return nullptr;
        }
        //Create stream
        IStream *stream;
        HRESULT hr = CreateStreamOnHGlobal(hmem,TRUE,&stream);
        if(FAILED(hr)){
            GlobalUnlock(hmem);
            GlobalFree(hmem);
            return nullptr;
        }
        return stream;
    }
    inline Bitmap LoadFromIStream(IStream *stream,int want){
        auto factory = GetFactory();
        if(factory == nullptr){
            return Bitmap();
        }
        ComPtr<IWICBitmapDecoder> decoder;
        HRESULT hr;
        hr = factory->CreateDecoderFromStream(
            stream,
            nullptr,
            WICDecodeMetadataCacheOnDemand,
            &decoder
        );
        if(FAILED(hr)){
            return Bitmap();
        }
        ComPtr<IWICBitmapFrameDecode> frame;
        hr = decoder->GetFrame(0, &frame);
        if(FAILED(hr)){
            return Bitmap();
        }
        ComPtr<IWICFormatConverter> converter;
        hr = factory->CreateFormatConverter(&converter);
        if(FAILED(hr)){
            return Bitmap();
        }
        const GUID *id;
        switch(want){
            case Format::Gray:
                id = &ID_Gray;
                break;
            case Format::GrayA:
                id = &ID_GrayA;
                break;
            case Format::RGB:
                id = &ID_RGB;
                break;
            case Format::RGBA:
                id = &ID_RGBA;
                break;
            default:
                return Bitmap();
        }
        hr = converter->Initialize(
            frame.Get(),
            *id,
            WICBitmapDitherTypeNone,
            nullptr,
            0.0,
            WICBitmapPaletteTypeMedianCut
        );
        if(FAILED(hr)){
            return Bitmap();
        }
        UINT w,h;
        hr = frame->GetSize(&w,&h);
        if(FAILED(hr)){
            return Bitmap();
        }
        //Alloc bitmap
        Bitmap bitmap;
        bitmap.w = w;
        bitmap.h = h;
        bitmap.channels = want;
        bitmap.pixels = AllocFrame(w,h,want);
        bitmap.pitch = w * want;
        //Copy pixels
        hr = converter->CopyPixels(
            nullptr,
            bitmap.pitch,
            bitmap.pitch * bitmap.h,
            (BYTE*)bitmap.pixels
        );
        if(FAILED(hr)){
            Free(bitmap.pixels);
            return Bitmap();
        }
        return bitmap;
    }
    inline Bitmap LoadFromCallback(const IOCallback &cb,int channals){
        ComPtr<IStream> stream(CallbackToIStream(cb));
        if(stream.Get() == nullptr){
            return Bitmap();
        }
        return LoadFromIStream(stream.Get(),channals);
    }
    inline Animation *LoadAnimation(const IOCallback &cb,int want){
        //TODO:
        ComPtr<IStream> stream(CallbackToIStream(cb));
        if(stream.Get() == nullptr){
            return nullptr;
        }
        auto factory = GetFactory();
        if(factory == nullptr){
            return nullptr;
        }
        ComPtr<IWICBitmapDecoder> decoder;
        HRESULT hr;
        hr = factory->CreateDecoderFromStream(
            stream.Get(),
            nullptr,
            WICDecodeMetadataCacheOnDemand,
            &decoder
        );
        if(FAILED(hr)){
            return nullptr;
        }
        UINT frame_count;
        UINT w,h;
        hr = decoder->GetFrameCount(&frame_count);
        if(FAILED(hr)){
            return nullptr;
        }
        ComPtr<IWICBitmapFrameDecode> frame;
        hr = decoder->GetFrame(0, &frame);
        if(FAILED(hr)){
            return nullptr;
        }
        hr = frame->GetSize(&w,&h);
        if(FAILED(hr)){
            return nullptr;
        }
        Animation *an = AllocAnimation(w,h,want,frame_count);
        UINT cur_frame = 0;

        do{
            ComPtr<IWICFormatConverter> converter;
            hr = factory->CreateFormatConverter(&converter);
            if(FAILED(hr)){
                goto err;
            }
            const GUID *id;
            switch(want){
                case Format::Gray:
                    id = &ID_Gray;
                    break;
                case Format::GrayA:
                    id = &ID_GrayA;
                    break;
                case Format::RGB:
                    id = &ID_RGB;
                    break;
                case Format::RGBA:
                    id = &ID_RGBA;
                    break;
                default:
                    goto err;
            }
            hr = converter->Initialize(
                frame.Get(),
                *id,
                WICBitmapDitherTypeNone,
                nullptr,
                0.0,
                WICBitmapPaletteTypeMedianCut
            );
            if(FAILED(hr)){
                goto err;
            }
            hr = converter->CopyPixels(
                nullptr,
                an->bitmaps[cur_frame].pitch,
                an->bitmaps[cur_frame].pitch * an->bitmaps[cur_frame].h,
                (BYTE*)an->bitmaps[cur_frame].pixels
            );
            if(FAILED(hr)){
                goto err;
            }
            ++cur_frame;
            if(cur_frame > frame_count){
                break;
            }
            //Get next frame
            hr = decoder->GetFrame(cur_frame, &frame);
            if(FAILED(hr)){
                goto err;
            }
        }
        while(true);

        return an;
        err:
            FreeAnimation(an);
            return nullptr;
    }

    #undef  _PIXUTIL_WIC_LOAD_ENTRY
    #define _PIXUTIL_WIC_LOAD_ENTRY Wincodec::LoadFromCallback,
} // namespace Wincodec
} // namespace PixImage

#endif // Wincodec

#if PIXUTIL_HAS_INCLUDE(<png.h>)
#include <png.h>
namespace PixImage{
namespace Png{
    inline Bitmap Load(const IOCallback &cb,int want){
        //Load from stream using libpng
        //TODO:
        return {};
    }
} // namespace Png
} // namespace PixImage

#endif // png

//Main function implementation
namespace PixImage{
    inline Bitmap LoadFromCallback(const IOCallback &cb,int want){
        using Loader = Bitmap(*)(const IOCallback &cb,int channals);
        Loader tables [] = {
            _PIXUTIL_WIC_LOAD_ENTRY
            nullptr
        };

        Bitmap bitmap;
        for(auto loader:tables){
            if(loader != nullptr){
                bitmap = loader(cb,want);
                if(!bitmap.empty()){
                    break;
                }
            }
        }
        return bitmap;
    }
    inline Bitmap LoadFromFilename(const char *filename,int want){
        FILE *fp;
        #ifdef _WIN32
        //Using unicode
        size_t len = strlen(filename);
        wchar_t *wfilename = (wchar_t *)PIXUTIL_MALLOC(sizeof(wchar_t) * (len + 1));
        MultiByteToWideChar(CP_UTF8,0,filename,-1,wfilename,len + 1);
        fp = _wfopen(wfilename,L"rb");
        PIXUTIL_FREE(wfilename);
        #else
        fp = fopen(filename,"rb");
        #endif
        if(fp == nullptr){
            return Bitmap();
        }
        auto cb = WrapFILE(fp);
        Bitmap bitmap = LoadFromCallback(cb,want);
        cb.close(cb.uptr);
        return bitmap;
    }

    inline IOCallback WrapFILE(FILE *f){
        IOCallback cb;
        cb.uptr = f;
        cb.read = [](void *uptr,void *ptr,size_t size) -> ssize_t{
            FILE *f = (FILE*)uptr;
            return fread(ptr,1,size,f);
        };
        cb.write = [](void *uptr,const void *ptr,size_t size) -> ssize_t{
            FILE *f = (FILE*)uptr;
            return fwrite(ptr,1,size,f);
        };
        cb.seek = [](void *uptr,int64_t offset,int whence) -> int64_t{
            FILE *f = (FILE*)uptr;
            return fseek(f,offset,whence);
        };
        cb.tell = [](void *uptr) -> int64_t{
            FILE *f = (FILE*)uptr;
            return ftell(f);
        };
        cb.close = [](void *uptr) -> void{
            FILE *f = (FILE*)uptr;
            fclose(f);
        };
        return cb;
    }
    inline IOCallback WrapMem(const void *mem,size_t n){
        auto cb = WrapMem(const_cast<void*>(mem),n);
        cb.write = [](void *uptr,const void *ptr,size_t size) -> ssize_t{
            return -1;
        };
        return cb;
    }
    inline IOCallback WrapMem(void *mem,size_t n){
        struct Recorder{
            void *mem;
            size_t n;
            size_t cur;
        };
        IOCallback cb;
        Recorder *rec = (Recorder*)Malloc(sizeof(Recorder));
        rec->mem = mem;
        rec->n = n;
        rec->cur = 0;
        cb.uptr = rec;
        cb.read = [](void *uptr,void *ptr,size_t size) -> ssize_t{
            Recorder *rec = (Recorder*)uptr;
            if(rec->cur + size > rec->n){
                size = rec->n - rec->cur;
            }
            memcpy(ptr,(char*)rec->mem + rec->cur,size);
            rec->cur += size;
            return size;
        };
        cb.write = [](void *uptr,const void *ptr,size_t size) -> ssize_t{
            Recorder *rec = (Recorder*)uptr;
            if(rec->cur + size > rec->n){
                size = rec->n - rec->cur;
            }
            memcpy((char*)rec->mem + rec->cur,ptr,size);
            rec->cur += size;
            return size;
        };
        cb.seek = [](void *uptr,int64_t offset,int whence) -> int64_t{
            Recorder *rec = (Recorder*)uptr;
            switch(whence){
                case SEEK_SET:
                    if(offset < 0 || offset > rec->n){
                        return -1;
                    }
                    rec->cur = offset;
                    break;
                case SEEK_CUR:
                    if(rec->cur + offset < 0 || rec->cur + offset > rec->n){
                        return -1;
                    }
                    rec->cur += offset;
                    break;
                case SEEK_END:
                    if(offset > 0 || offset + rec->n < 0){
                        return -1;
                    }
                    rec->cur = rec->n + offset;
                    break;
                default:
                    return -1;
            }
            return rec->cur;
        };
        cb.tell = [](void *uptr) -> int64_t{
            Recorder *rec = (Recorder*)uptr;
            return rec->cur;
        };
        cb.close = [](void *uptr) -> void{
            Recorder *rec = (Recorder*)uptr;
            Free(rec);
        };
        return cb;
    }
} // namespace PixImage

#endif // _BTK_PROJECT_PIXUTIL_IMAGE_HPP_
