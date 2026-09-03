#pragma once

// Replaces legacy/parameters.cpp's loose globals (`double U=3, mu=0., ...;
// int Lx=8, Ly=8, N=Lx*Ly; ...`). Defaults below match that file exactly.
//
// N, Ntau4 and Ntau were separately-stored globals in the original
// (`N=Lx*Ly`, `Ntau4=10*Nw4`, `Ntau=4*Nw`) with nothing preventing Lx/Ly or
// Nw4/Nw from being changed afterwards and leaving them stale. Here they're
// computed from the fields they derive from, so that can't happen.
namespace dmft {

struct Parameters {
    double U = 3.0;
    double mu = 0.0;
    double t = 1.0;
    double t_1 = 0.0;
    double t_2 = 0.0;
    double beta = 4.0;

    int Lx = 8;
    int Ly = 8;
    int N() const { return Lx * Ly; }

    int Nbath_raw = 3;
    int Nbath_fine = 9;
    double Bath_fine_halfwidth = 3.0;

    int NwED = 20;
    int Nw = 200;
    int Nw4 = 16;
    int Ntau4() const { return 10 * Nw4; }
    int Ntau() const { return 4 * Nw; }
};

}  // namespace dmft
