# THIS FILE HAS BEEN AUTOMATICALLY GENERATED.

solution_dict = {

    "domain_extents" : (1.0, 1.2, 0.8),

    "pressure1_scale" : 1.0,

    "saturation2_scale" : 0.01,

    "GRAVITY_MAGNITUDE" : "9.8",

    "VISCOSITY1" : "1.725e-05",

    "VISCOSITY2" : "0.001",

    "DENSITY1" : "1.284",

    "DENSITY2" : "1000.0",

    "POROSITY" : "0.4",

    "ABSOLUTE_PERMEABILITY" : "1.567346939e-09",

    "GRAVITY_DIRECTION_1D" : "-1.",

    "PRESSURE_1D" : "1.0*sin(pi*X[0])**2",

    "SATURATION2_1D" : "0.0225*X[0]**2*(-3.0*X[0] + 3.0)",

    "DARCY_VELOCITY1_MAGNITUDE_1D" : "9.0860692115942e-5*sqrt((-0.0225*X[0]**2*(-3.0*X[0] + 3.0) + 1.)**4*(2.0*pi*sin(pi*X[0])*cos(pi*X[0]) + 12.5832)**2)",

    "SOURCE_SATURATION1_1D" : "-9.0860692115942e-5*(0.135*X[0]**2 - 0.09*X[0]*(-3.0*X[0] + 3.0))*(-0.0225*X[0]**2*(-3.0*X[0] + 3.0) + 1.)*(2.0*pi*sin(pi*X[0])*cos(pi*X[0]) + 12.5832) - 9.0860692115942e-5*(-2.0*pi**2*sin(pi*X[0])**2 + 2.0*pi**2*cos(pi*X[0])**2)*(-0.0225*X[0]**2*(-3.0*X[0] + 3.0) + 1.)**2",

    "DARCY_VELOCITY2_MAGNITUDE_1D" : "7.9346938786875e-10*sqrt(X[0]**8*(-3.0*X[0] + 3.0)**4*(2.0*pi*sin(pi*X[0])*cos(pi*X[0]) + 9800.0)**2)",

    "SOURCE_SATURATION2_1D" : "-7.9346938786875e-10*X[0]**4*(-3.0*X[0] + 3.0)**2*(-2.0*pi**2*sin(pi*X[0])**2 + 2.0*pi**2*cos(pi*X[0])**2) - 7.9346938786875e-10*X[0]**4*(18.0*X[0] - 18.0)*(2.0*pi*sin(pi*X[0])*cos(pi*X[0]) + 9800.0) - 3.173877551475e-9*X[0]**3*(-3.0*X[0] + 3.0)**2*(2.0*pi*sin(pi*X[0])*cos(pi*X[0]) + 9800.0)",

    "GRAVITY_DIRECTION_2D" : "-1. 0. ",

    "PRESSURE_2D" : "1.0*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2",

    "SATURATION2_2D" : "0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)",

    "DARCY_VELOCITY1_MAGNITUDE_2D" : "sqrt(8.25566537178801e-9*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*cos(pi*X[0]) + 12.5832)**2*(-0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) + 1.)**4 + 2.29324038105222e-8*pi**2*(-0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) + 1.)**4*sin(pi*X[0])**4*sin(0.833333333333333*pi*X[1])**2*cos(0.833333333333333*pi*X[1])**2)",

    "SOURCE_SATURATION1_2D" : "-9.0860692115942e-5*(-2.0*pi**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2 + 2.0*pi**2*sin(0.833333333333333*pi*X[1])**2*cos(pi*X[0])**2)*(-0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) + 1.)**2 - 0.000151434486859903*pi*(0.17578125*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0) - 0.140625*X[0]**2*X[1]*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0))*(-0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) + 1.)*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])*cos(0.833333333333333*pi*X[1]) - 9.0860692115942e-5*(0.2109375*X[0]**2*X[1]**2*(-2.5*X[1] + 3.0) - 0.140625*X[0]*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0))*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*cos(pi*X[0]) + 12.5832)*(-0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) + 1.) + 0.000126195405716586*pi**2*(-0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) + 1.)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2 - 0.000126195405716586*pi**2*(-0.03515625*X[0]**2*X[1]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) + 1.)**2*sin(pi*X[0])**2*cos(0.833333333333333*pi*X[1])**2",

    "DARCY_VELOCITY2_MAGNITUDE_2D" : "sqrt(3.75267070224767e-18*X[0]**8*X[1]**8*(-3.0*X[0] + 3.0)**4*(-2.5*X[1] + 3.0)**4*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*cos(pi*X[0]) + 9800.0)**2 + 1.04240852840213e-17*pi**2*X[0]**8*X[1]**8*(-3.0*X[0] + 3.0)**4*(-2.5*X[1] + 3.0)**4*sin(pi*X[0])**4*sin(0.833333333333333*pi*X[1])**2*cos(0.833333333333333*pi*X[1])**2)",

    "SOURCE_SATURATION2_2D" : "-1.93718112272644e-9*X[0]**4*X[1]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-2.0*pi**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2 + 2.0*pi**2*sin(0.833333333333333*pi*X[1])**2*cos(pi*X[0])**2) + 2.69052933712006e-9*pi**2*X[0]**4*X[1]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2 - 2.69052933712006e-9*pi**2*X[0]**4*X[1]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*sin(pi*X[0])**2*cos(0.833333333333333*pi*X[1])**2 - 3.22863520454407e-9*pi*X[0]**4*X[1]**4*(-3.0*X[0] + 3.0)**2*(12.5*X[1] - 15.0)*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])*cos(0.833333333333333*pi*X[1]) - 1.93718112272644e-9*X[0]**4*X[1]**4*(18.0*X[0] - 18.0)*(-2.5*X[1] + 3.0)**2*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*cos(pi*X[0]) + 9800.0) - 1.29145408181763e-8*pi*X[0]**4*X[1]**3*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])*cos(0.833333333333333*pi*X[1]) - 7.74872449090576e-9*X[0]**3*X[1]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*cos(pi*X[0]) + 9800.0)",

    "GRAVITY_DIRECTION_3D" : "-1. 0. 0. ",

    "PRESSURE_3D" : "1.0*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2",

    "SATURATION2_3D" : "0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0)",

    "DARCY_VELOCITY1_MAGNITUDE_3D" : "sqrt(8.25566537178801e-9*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2*cos(pi*X[0]) + 12.5832)**2*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)**4 + 5.15979085736751e-8*pi**2*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)**4*sin(pi*X[0])**4*sin(0.833333333333333*pi*X[1])**4*sin(1.25*pi*X[2])**2*cos(1.25*pi*X[2])**2 + 2.29324038105222e-8*pi**2*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)**4*sin(pi*X[0])**4*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**4*cos(0.833333333333333*pi*X[1])**2)",

    "SOURCE_SATURATION1_3D" : "-9.0860692115942e-5*(-2.0*pi**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2 + 2.0*pi**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2*cos(pi*X[0])**2)*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)**2 - 9.0860692115942e-5*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2*cos(pi*X[0]) + 12.5832)*(0.7415771484375*X[0]**2*X[1]**2*X[2]**2*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) - 0.494384765625*X[0]*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0))*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.) - 0.000227151730289855*pi*(0.926971435546875*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0) - 0.494384765625*X[0]**2*X[1]**2*X[2]*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0))*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])*cos(1.25*pi*X[2]) - 0.000151434486859903*pi*(0.61798095703125*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-3.75*X[2] + 3.0) - 0.494384765625*X[0]**2*X[1]*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0))*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])*sin(1.25*pi*X[2])**2*cos(0.833333333333333*pi*X[1]) + 0.000410135068578905*pi**2*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2 - 0.000283939662862319*pi**2*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*cos(1.25*pi*X[2])**2 - 0.000126195405716586*pi**2*(-0.12359619140625*X[0]**2*X[1]**2*X[2]**2*(-3.0*X[0] + 3.0)*(-2.5*X[1] + 3.0)*(-3.75*X[2] + 3.0) + 1.)**2*sin(pi*X[0])**2*sin(1.25*pi*X[2])**2*cos(0.833333333333333*pi*X[1])**2",

    "DARCY_VELOCITY2_MAGNITUDE_3D" : "sqrt(5.73258671850129e-16*X[0]**8*X[1]**8*X[2]**8*(-3.0*X[0] + 3.0)**4*(-2.5*X[1] + 3.0)**4*(-3.75*X[2] + 3.0)**4*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2*cos(pi*X[0]) + 9800.0)**2 + 3.5828666990633e-15*pi**2*X[0]**8*X[1]**8*X[2]**8*(-3.0*X[0] + 3.0)**4*(-2.5*X[1] + 3.0)**4*(-3.75*X[2] + 3.0)**4*sin(pi*X[0])**4*sin(0.833333333333333*pi*X[1])**4*sin(1.25*pi*X[2])**2*cos(1.25*pi*X[2])**2 + 1.59238519958369e-15*pi**2*X[0]**8*X[1]**8*X[2]**8*(-3.0*X[0] + 3.0)**4*(-2.5*X[1] + 3.0)**4*(-3.75*X[2] + 3.0)**4*sin(pi*X[0])**4*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**4*cos(0.833333333333333*pi*X[1])**2)",

    "SOURCE_SATURATION2_3D" : "-2.39428208833071e-8*X[0]**4*X[1]**4*X[2]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*(-2.0*pi**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2 + 2.0*pi**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2*cos(pi*X[0])**2) + 1.08075233153817e-7*pi**2*X[0]**4*X[1]**4*X[2]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2 - 7.48213152603348e-8*pi**2*X[0]**4*X[1]**4*X[2]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*cos(1.25*pi*X[2])**2 - 3.32539178934821e-8*pi**2*X[0]**4*X[1]**4*X[2]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*sin(pi*X[0])**2*sin(1.25*pi*X[2])**2*cos(0.833333333333333*pi*X[1])**2 - 5.98570522082678e-8*pi*X[0]**4*X[1]**4*X[2]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(28.125*X[2] - 22.5)*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])*cos(1.25*pi*X[2]) - 3.99047014721786e-8*pi*X[0]**4*X[1]**4*X[2]**4*(-3.0*X[0] + 3.0)**2*(12.5*X[1] - 15.0)*(-3.75*X[2] + 3.0)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])*sin(1.25*pi*X[2])**2*cos(0.833333333333333*pi*X[1]) - 2.39428208833071e-8*X[0]**4*X[1]**4*X[2]**4*(18.0*X[0] - 18.0)*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2*cos(pi*X[0]) + 9800.0) - 2.39428208833071e-7*pi*X[0]**4*X[1]**4*X[2]**3*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])*cos(1.25*pi*X[2]) - 1.59618805888714e-7*pi*X[0]**4*X[1]**3*X[2]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*sin(pi*X[0])**2*sin(0.833333333333333*pi*X[1])*sin(1.25*pi*X[2])**2*cos(0.833333333333333*pi*X[1]) - 9.57712835332286e-8*X[0]**3*X[1]**4*X[2]**4*(-3.0*X[0] + 3.0)**2*(-2.5*X[1] + 3.0)**2*(-3.75*X[2] + 3.0)**2*(2.0*pi*sin(pi*X[0])*sin(0.833333333333333*pi*X[1])**2*sin(1.25*pi*X[2])**2*cos(pi*X[0]) + 9800.0)",

    }
