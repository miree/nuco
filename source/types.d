import nucd.geometry;

// some types
alias Vector!(2,real) Vec2;
alias Vector!(3,real) Vec3;
alias Matrix!(3,3,rowMajor,real) Mat3;



// some constants
immutable real c   = 299.792458; // speed of light   [fm/zs]
immutable real ahc = 1.43996442; // alpha*hbar*c     [MeV*fm]
immutable real u   = 931.49406;  // atomic mass unit [MeV/c^2]
immutable real mp  = 938.272013; // proton  mass     [MeV/c^2]
immutable real mn  = 939.565346; // neutron mass     [MeV/c^2]
immutable real me  = 0.51099893; // electron mass    [MeV/c^2]
