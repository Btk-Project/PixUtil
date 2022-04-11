#define PIXUTIL_SDL_EXTERNAL
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <iostream>
#include <ctime>
#define PIXUTIL_MATRIX_IO
#include "../pixutil_filter.hpp"
#include "../pixutil.hpp"

int main(int argc,char **argv){
    SDL_Surface *surf = IMG_Load(
        argv[1]
    );
    if(surf == nullptr){
        printf("Load failed\n");
        return -1;
    }
    SDL_Surface *out = SDL_DuplicateSurface(surf);
    SDL_Surface *out1 = SDL_DuplicateSurface(surf);
    //Gray buffer
    PixUtil::SDLSurfaceView view(surf);
    PixUtil::SDLSurfaceView dst(out);
    PixUtil::SDLSurfaceView dst1(out1);

    PixFilter::FFTMatrix mat[3];
    clock_t prev = clock();
    PixFilter::FFT2D(mat,view);
    clock_t cur = clock();
    std::cout << "FFT2D: " << (cur - prev) / CLOCKS_PER_SEC << "s" << std::endl;

    prev = clock();
    PixFilter::IFFT2D(dst,mat);
    cur = clock();
    std::cout << "IFFT2D: " << (cur - prev) / CLOCKS_PER_SEC << "s" << std::endl;


    PixFilter::MergeFFTMatrixToView(dst1,mat,3);

    PixUtil::SDLShowView("IFFT",dst);
    PixUtil::SDLShowView("FFT",dst1);
    PixUtil::SDLShowView("Org",view);

    PixUtil::SDLWaitForQuit();

    SDL_FreeSurface(surf);
    SDL_FreeSurface(out);
    SDL_FreeSurface(out1);
    SDL_Quit();
}