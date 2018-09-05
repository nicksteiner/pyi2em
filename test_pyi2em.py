import pyi2em

def main()
    print("FREQ: 3.0 [GHZ], CDC: 11.3 + i1.5, CL: 10 [cm], RMS: .25 [cm], INC: 30 [deg], CORRF: Gaussian")
    e_ = pyi2em.test_I2EM()
    print("Emissivity: %g [V], %g [H]".formaty(*e_))

if __name__ == '__main__':
    main()