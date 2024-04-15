module GPL(input A, B, output P, G);

    xor (P, A, B);
    and (G, A, B);

endmodule

module GrayCell(input Gik, Pik, Gk1j, output Gij);

    wire net1;
    and (net1, Pik, Gk1j);
    or (Gij, net1, Gik);

endmodule

module BlackCell(input Gik, Pik, Gk1j, Pk1j, output Gij, Pij);

    wire net1;
    and (net1, Pik, Gk1j);
    or (Gij, net1, Gik);
    and (Pij, Pik, Pk1j);

endmodule

module KSA16(input [15:0] A, B, input Cin, output [15:0] S, output Cout);

    wire [15:0] G;//From individual PG logic
    wire [15:0] P;
    wire [13:0] L0G;//Coming out of level 0 group PG G bus
    wire [13:0] L0P;//Coming out of level 0 group PG P bus
    wire [11:0] L1G;//Coming out of level 1 group PG G bus
    wire [11:0] L1P;//Coming out of level 1 group PG P bus
    wire [7:0] L2G;//Coming out of level 2 group PG G bus
    wire [7:0] L2P;//Coming out of level 2 group PG P bus
    wire [14:0] G_i0;//To Sum logic, coming out of level 3 group PG logic

    //Individual PG Logic
    GPL gp0(.A(A[0]),.B(B[0]),.P(P[0]),.G(G[0])); 
    GPL gp1(.A(A[1]),.B(B[1]),.P(P[1]),.G(G[1])); 
    GPL gp2(.A(A[2]),.B(B[2]),.P(P[2]),.G(G[2]));
    GPL gp3(.A(A[3]),.B(B[3]),.P(P[3]),.G(G[3]));
    GPL gp4(.A(A[4]),.B(B[4]),.P(P[4]),.G(G[4]));
    GPL gp5(.A(A[5]),.B(B[5]),.P(P[5]),.G(G[5]));
    GPL gp6(.A(A[6]),.B(B[6]),.P(P[6]),.G(G[6]));
    GPL gp7(.A(A[7]),.B(B[7]),.P(P[7]),.G(G[7]));
    GPL gp8(.A(A[8]),.B(B[8]),.P(P[8]),.G(G[8])); 
    GPL gp9(.A(A[9]),.B(B[9]),.P(P[9]),.G(G[9])); 
    GPL gp10(.A(A[10]),.B(B[10]),.P(P[10]),.G(G[10]));
    GPL gp11(.A(A[11]),.B(B[11]),.P(P[11]),.G(G[11]));
    GPL gp12(.A(A[12]),.B(B[12]),.P(P[12]),.G(G[12]));
    GPL gp13(.A(A[13]),.B(B[13]),.P(P[13]),.G(G[13]));
    GPL gp14(.A(A[14]),.B(B[14]),.P(P[14]),.G(G[14]));
    GPL gp15(.A(A[15]),.B(B[15]),.P(P[15]),.G(G[15]));  

    //Level 0 Group PG Logic
    GrayCell gray01(.Gik(G[0]),.Pik(P[0]),.Gk1j(Cin),.Gij(G_i0[0]));
    BlackCell black02(.Gik(G[1]),.Pik(P[1]),.Gk1j(G[0]),.Pk1j(P[0]),.Gij(L0G[0]),.Pij(L0P[0]));
    BlackCell black03(.Gik(G[2]),.Pik(P[2]),.Gk1j(G[1]),.Pk1j(P[1]),.Gij(L0G[1]),.Pij(L0P[1]));
    BlackCell black04(.Gik(G[3]),.Pik(P[3]),.Gk1j(G[2]),.Pk1j(P[2]),.Gij(L0G[2]),.Pij(L0P[2]));
    BlackCell black05(.Gik(G[4]),.Pik(P[4]),.Gk1j(G[3]),.Pk1j(P[3]),.Gij(L0G[3]),.Pij(L0P[3]));
    BlackCell black06(.Gik(G[5]),.Pik(P[5]),.Gk1j(G[4]),.Pk1j(P[4]),.Gij(L0G[4]),.Pij(L0P[4]));
    BlackCell black07(.Gik(G[6]),.Pik(P[6]),.Gk1j(G[5]),.Pk1j(P[5]),.Gij(L0G[5]),.Pij(L0P[5]));
    BlackCell black08(.Gik(G[7]),.Pik(P[7]),.Gk1j(G[6]),.Pk1j(P[6]),.Gij(L0G[6]),.Pij(L0P[6]));
    BlackCell black09(.Gik(G[8]),.Pik(P[8]),.Gk1j(G[7]),.Pk1j(P[7]),.Gij(L0G[7]),.Pij(L0P[7]));
    BlackCell black010(.Gik(G[9]),.Pik(P[9]),.Gk1j(G[8]),.Pk1j(P[8]),.Gij(L0G[8]),.Pij(L0P[8]));
    BlackCell black011(.Gik(G[10]),.Pik(P[10]),.Gk1j(G[9]),.Pk1j(P[9]),.Gij(L0G[9]),.Pij(L0P[9]));
    BlackCell black012(.Gik(G[11]),.Pik(P[11]),.Gk1j(G[10]),.Pk1j(P[10]),.Gij(L0G[10]),.Pij(L0P[10]));
    BlackCell black013(.Gik(G[12]),.Pik(P[12]),.Gk1j(G[11]),.Pk1j(P[11]),.Gij(L0G[11]),.Pij(L0P[11]));
    BlackCell black014(.Gik(G[13]),.Pik(P[13]),.Gk1j(G[12]),.Pk1j(P[12]),.Gij(L0G[12]),.Pij(L0P[12]));
    BlackCell black015(.Gik(G[14]),.Pik(P[14]),.Gk1j(G[13]),.Pk1j(P[13]),.Gij(L0G[13]),.Pij(L0P[13]));

    //Level 1 Group PG Logic
    GrayCell gray12(.Gik(L0G[0]),.Pik(L0P[0]),.Gk1j(Cin),.Gij(G_i0[1]));
    GrayCell gray13(.Gik(L0G[1]),.Pik(L0P[1]),.Gk1j(G_i0[0]),.Gij(G_i0[2]));
    BlackCell black14(.Gik(L0G[2]),.Pik(L0P[2]),.Gk1j(L0G[0]),.Pk1j(L0P[0]),.Gij(L1G[0]),.Pij(L1P[0]));
    BlackCell black15(.Gik(L0G[3]),.Pik(L0P[3]),.Gk1j(L0G[1]),.Pk1j(L0P[1]),.Gij(L1G[1]),.Pij(L1P[1]));
    BlackCell black16(.Gik(L0G[4]),.Pik(L0P[4]),.Gk1j(L0G[2]),.Pk1j(L0P[2]),.Gij(L1G[2]),.Pij(L1P[2]));
    BlackCell black17(.Gik(L0G[5]),.Pik(L0P[5]),.Gk1j(L0G[3]),.Pk1j(L0P[3]),.Gij(L1G[3]),.Pij(L1P[3]));
    BlackCell black18(.Gik(L0G[6]),.Pik(L0P[6]),.Gk1j(L0G[4]),.Pk1j(L0P[4]),.Gij(L1G[4]),.Pij(L1P[4]));
    BlackCell black19(.Gik(L0G[7]),.Pik(L0P[7]),.Gk1j(L0G[5]),.Pk1j(L0P[5]),.Gij(L1G[5]),.Pij(L1P[5]));
    BlackCell black110(.Gik(L0G[8]),.Pik(L0P[8]),.Gk1j(L0G[6]),.Pk1j(L0P[6]),.Gij(L1G[6]),.Pij(L1P[6]));
    BlackCell black111(.Gik(L0G[9]),.Pik(L0P[9]),.Gk1j(L0G[7]),.Pk1j(L0P[7]),.Gij(L1G[7]),.Pij(L1P[7]));
    BlackCell black112(.Gik(L0G[10]),.Pik(L0P[10]),.Gk1j(L0G[8]),.Pk1j(L0P[8]),.Gij(L1G[8]),.Pij(L1P[8]));
    BlackCell black113(.Gik(L0G[11]),.Pik(L0P[11]),.Gk1j(L0G[9]),.Pk1j(L0P[9]),.Gij(L1G[9]),.Pij(L1P[9]));
    BlackCell black114(.Gik(L0G[12]),.Pik(L0P[12]),.Gk1j(L0G[10]),.Pk1j(L0P[10]),.Gij(L1G[10]),.Pij(L1P[10]));
    BlackCell black115(.Gik(L0G[13]),.Pik(L0P[13]),.Gk1j(L0G[11]),.Pk1j(L0P[11]),.Gij(L1G[11]),.Pij(L1P[11]));

    //Level 2 Group PG Logic
    GrayCell gray24(.Gik(L1G[0]),.Pik(L1P[0]),.Gk1j(Cin),.Gij(G_i0[3]));
    GrayCell gray25(.Gik(L1G[1]),.Pik(L1P[1]),.Gk1j(G_i0[0]),.Gij(G_i0[4]));
    GrayCell gray26(.Gik(L1G[2]),.Pik(L1P[2]),.Gk1j(G_i0[1]),.Gij(G_i0[5]));
    GrayCell gray27(.Gik(L1G[3]),.Pik(L1P[3]),.Gk1j(G_i0[2]),.Gij(G_i0[6]));
    BlackCell black28(.Gik(L1G[4]),.Pik(L1P[4]),.Gk1j(L1G[0]),.Pk1j(L1P[0]),.Gij(L2G[0]),.Pij(L2P[0]));
    BlackCell black29(.Gik(L1G[5]),.Pik(L1P[5]),.Gk1j(L1G[1]),.Pk1j(L1P[1]),.Gij(L2G[1]),.Pij(L2P[1]));
    BlackCell black210(.Gik(L1G[6]),.Pik(L1P[6]),.Gk1j(L1G[2]),.Pk1j(L1P[2]),.Gij(L2G[2]),.Pij(L2P[2]));
    BlackCell black211(.Gik(L1G[7]),.Pik(L1P[7]),.Gk1j(L1G[3]),.Pk1j(L1P[3]),.Gij(L2G[3]),.Pij(L2P[3]));
    BlackCell black212(.Gik(L1G[8]),.Pik(L1P[8]),.Gk1j(L1G[4]),.Pk1j(L1P[4]),.Gij(L2G[4]),.Pij(L2P[4]));
    BlackCell black213(.Gik(L1G[9]),.Pik(L1P[9]),.Gk1j(L1G[5]),.Pk1j(L1P[5]),.Gij(L2G[5]),.Pij(L2P[5]));
    BlackCell black214(.Gik(L1G[10]),.Pik(L1P[10]),.Gk1j(L1G[6]),.Pk1j(L1P[6]),.Gij(L2G[6]),.Pij(L2P[6]));
    BlackCell black215(.Gik(L1G[11]),.Pik(L1P[11]),.Gk1j(L1G[7]),.Pk1j(L1P[7]),.Gij(L2G[7]),.Pij(L2P[7]));

    //Level 3 Group PG Logic
    GrayCell gray38(.Gik(L2G[0]),.Pik(L2P[0]),.Gk1j(Cin),.Gij(G_i0[7]));
    GrayCell gray39(.Gik(L2G[1]),.Pik(L2P[1]),.Gk1j(G_i0[0]),.Gij(G_i0[8]));
    GrayCell gray310(.Gik(L2G[2]),.Pik(L2P[2]),.Gk1j(G_i0[1]),.Gij(G_i0[9]));
    GrayCell gray311(.Gik(L2G[3]),.Pik(L2P[3]),.Gk1j(G_i0[3]),.Gij(G_i0[10]));
    GrayCell gray312(.Gik(L2G[4]),.Pik(L2P[4]),.Gk1j(G_i0[4]),.Gij(G_i0[11]));
    GrayCell gray313(.Gik(L2G[5]),.Pik(L2P[5]),.Gk1j(G_i0[5]),.Gij(G_i0[12]));
    GrayCell gray314(.Gik(L2G[6]),.Pik(L2P[6]),.Gk1j(G_i0[6]),.Gij(G_i0[13]));
    GrayCell gray315(.Gik(L2G[7]),.Pik(L2P[7]),.Gk1j(G_i0[7]),.Gij(G_i0[14]));

    //Sum Logic
    xor (S[0], Cin, P[0]);
    xor (S[1], G_i0[0], P[1]);
    xor (S[2], G_i0[1], P[2]);
    xor (S[3], G_i0[2], P[3]);
    xor (S[4], G_i0[3], P[4]);
    xor (S[5], G_i0[4], P[5]);
    xor (S[6], G_i0[5], P[6]);
    xor (S[7], G_i0[6], P[7]);
    xor (S[8], G_i0[7], P[8]);
    xor (S[9], G_i0[8], P[9]);
    xor (S[10], G_i0[9], P[10]);
    xor (S[11], G_i0[10], P[11]);
    xor (S[12], G_i0[11], P[12]);
    xor (S[13], G_i0[12], P[13]);
    xor (S[14], G_i0[13], P[14]);
    xor (S[15], G_i0[14], P[15]);
    GrayCell grayOUT(.Gik(G[15]), .Pik(P[15]), .Gk1j(C_in), .Gij(Cout));

endmodule