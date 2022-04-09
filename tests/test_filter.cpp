#include <SDL2/SDL.h>

#include <iostream>

#define PIXUTIL_MATRIX_IO
#define PIXUTIL_SDL_EXTERNAL

#include "../pixutil_filter.hpp"
#include "../pixutil_image.hpp"

using namespace PixFilter;
using namespace PixUtil;

#if __has_include(<SDL2/SDL_image.h>)
    #include <SDL2/SDL_image.h>
#endif

int main(int argc,char **argv){
    int w = 1000;
    int h = 600;
    if(argc == 2){
        //WINDOW w h
    }

    auto mat1 = EdgeDetectFilter();

    std::cout << "mat1: " << mat1 << std::endl;


    //Print mat1
    SDL_Window *win;
    SDL_Renderer *render;
    SDL_Texture *tex = nullptr;
    SDL_CreateWindowAndRenderer(w,h,SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE,&win,&render);

    SDL_Event event;
    while(SDL_WaitEvent(&event)){
        if(event.type == SDL_QUIT){
            break;
        }
        if(event.type == SDL_DROPFILE){
            //Drop test file
            #if __has_include(<SDL2/SDL_image.h>)
            SDL_Surface *surf = IMG_Load(event.drop.file);
            #else
            auto bitmap = PixImage::LoadFromFilename(event.drop.file);

            if(bitmap.empty()){
                continue;
            }
            SDL_Surface *surf;
            surf = SDL_CreateRGBSurfaceWithFormatFrom(
                bitmap.pixels,
                bitmap.w,
                bitmap.h,
                32,
                bitmap.pitch,
                SDL_PIXELFORMAT_RGBA32
            );
            if(surf == nullptr){
                printf("Load failed\n");
                continue;
            }
            #endif
            if(SDL_ISPIXELFORMAT_ALPHA(surf->format->format)){
                //Alpha
                printf("Alpha format\n");
            }
            printf("format is %s\n",SDL_GetPixelFormatName(surf->format->format));
            printf("surface pitch is %d\n",surf->pitch);
            printf("surface w is %d\n",surf->w);
            printf("surface w is %d\n",surf->h);


            SDL_Surface *dst = SDL_DuplicateSurface(surf);
            SDLSurfaceView view(surf);
            SDLSurfaceView dst_view(dst);

            //Test locating
            {
                SDL_Surface *new_surf = SDL_ConvertSurfaceFormat(
                    surf,
                    SDL_PIXELFORMAT_RGBA32,
                    0
                );
                SDLSurfaceView new_view(new_surf);

                //Check color is same 
                for(int y = 0;y < surf->h;y++){
                    for(int x = 0;x < surf->w;x++){
                        Color c1 = view.at(x,y);
                        Color c2 = new_view.at(x,y);
                        assert(c1.r == c2.r);
                        assert(c1.g == c2.g);
                        assert(c1.b == c2.b);
                        assert(c1.a == c2.a);
                    }
                }

                SDL_FreeSurface(new_surf);
            }

            printf("Begin filter\n");
            Filter2D(dst_view,view,mat1);
            printf("End filter\n");

            SDL_DestroyTexture(tex);
            tex = SDL_CreateTextureFromSurface(render,dst);
            SDL_FreeSurface(surf);
            SDL_FreeSurface(dst);

            #if !__has_include(<SDL2/SDL_image.h>)
            PixImage::Free(bitmap.pixels);
            #endif
        }

        SDL_RenderClear(render);

        if(tex != nullptr){
            int win_w,win_h;
            int tex_w,tex_h;

            SDL_QueryTexture(tex,nullptr,nullptr,&tex_w,&tex_h);
            SDL_GetWindowSize(win,&win_w,&win_h);
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

            SDL_RenderCopy(render,tex,nullptr,&rect);
        }
        SDL_RenderPresent(render);
    }

    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(render);
    SDL_DestroyWindow(win);

    SDL_Quit();
}