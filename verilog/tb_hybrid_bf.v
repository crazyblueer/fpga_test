`timescale 1ns/1ps
module tb_hybrid_bf;
    parameter M = 200;
    reg clk = 0, rst = 1, start = 0;
    always #10 clk = ~clk;

    reg x_we = 0; reg [5:0] x_addr = 0; reg signed [15:0] x_re = 0, x_im = 0;
    reg c_we = 0; reg [5:0] c_addr = 0; reg signed [15:0] c_re = 0, c_im = 0;
    reg b_we = 0; reg [2:0] b_addr = 0; reg signed [15:0] b_re = 0, b_im = 0;
    wire signed [23:0] y_re, y_im; wire done;

    hybrid_bf dut (clk, rst, start,
        x_we, x_addr, x_re, x_im,
        c_we, c_addr, c_re, c_im,
        b_we, b_addr, b_re, b_im,
        y_re, y_im, done);

    reg [31:0] xmem [0:64*M-1];
    reg [31:0] cmem [0:63];
    reg [31:0] bmem [0:7];
    reg [47:0] gmem [0:M-1];
    integer t, k, errs;

    initial begin
        $readmemh("x_vec.hex", xmem);  $readmemh("coef.hex", cmem);
        $readmemh("bb.hex", bmem);     $readmemh("y_gold.hex", gmem);
        errs = 0;
        repeat (3) @(posedge clk); rst = 0;

        for (k = 0; k < 64; k = k + 1) begin
            @(posedge clk); c_we <= 1; c_addr <= k; {c_re, c_im} <= cmem[k];
        end
        for (k = 0; k < 8; k = k + 1) begin
            @(posedge clk); c_we <= 0; b_we <= 1; b_addr <= k; {b_re, b_im} <= bmem[k];
        end
        @(posedge clk); b_we <= 0;

        for (t = 0; t < M; t = t + 1) begin
            for (k = 0; k < 64; k = k + 1) begin
                @(posedge clk); x_we <= 1; x_addr <= k; {x_re, x_im} <= xmem[t*64 + k];
            end
            @(posedge clk); x_we <= 0; start <= 1;
            @(posedge clk); start <= 0;
            @(posedge clk);
            wait (done == 1);
            @(posedge clk);
            if ({y_re, y_im} !== gmem[t]) begin
                errs = errs + 1;
                if (errs < 10)
                    $display("MISMATCH t=%0d  hw=%h  gold=%h", t, {y_re, y_im}, gmem[t]);
            end
        end
        $display("DONE: %0d snapshots, %0d errors", M, errs);
        $finish;
    end
endmodule