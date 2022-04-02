#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#define PIXUTIL_SDL_EXTERNAL

#include "../pixutil_filter.hpp"

using namespace PixFilter;
using namespace PixUtil;

int main(int argc,char **argv){
    int w = 1000;
    int h = 600;
    if(argc == 2){
        //WINDOW w h
    }

    auto mat1 = ImportMatrix<double,3,3>(1, 1, 1, 1, -8, 1, 1, 1, 1);

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
            SDL_Surface *surf = IMG_Load(event.drop.file);
            if(surf == nullptr){
                printf("Load failed\n");
                continue;
            }
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
        }

        SDL_RenderClear(render);

        if(tex != nullptr){
            SDL_RenderCopy(render,tex,nullptr,nullptr);
        }
        SDL_RenderPresent(render);
    }

    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(render);
    SDL_DestroyWindow(win);

    SDL_Quit();
}