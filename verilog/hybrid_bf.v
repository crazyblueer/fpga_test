module hybrid_bf (
    input  clk, rst, start,
    // write ports (later driven by the ARM)
    input  x_we, input [5:0] x_addr, input signed [15:0] x_re, x_im,
    input  c_we, input [5:0] c_addr, input signed [15:0] c_re, c_im,
    input  b_we, input [2:0] b_addr, input signed [15:0] b_re, b_im,
    output reg signed [23:0] y_re, y_im,
    output reg done
);
    localparam IDLE = 2'd0, S1 = 2'd1, S2 = 2'd2, FIN = 2'd3;
    reg [1:0] st;
    reg [6:0] cnt;

    reg signed [15:0] xr [0:63], xi [0:63];
    reg signed [15:0] cr [0:63], ci [0:63];
    reg signed [15:0] br [0:7],  bi [0:7];
    reg signed [19:0] zr [0:7],  zi [0:7];     // analog-stage outputs (Q4.15)
    reg signed [39:0] accr, acci;

    always @(posedge clk) begin
        if (x_we) begin xr[x_addr] <= x_re; xi[x_addr] <= x_im; end
        if (c_we) begin cr[c_addr] <= c_re; ci[c_addr] <= c_im; end
        if (b_we) begin br[b_addr] <= b_re; bi[b_addr] <= b_im; end
    end

    wire [5:0] i = cnt[5:0];
    wire [2:0] j = cnt[2:0];

    // stage 1: analog (phase-only) combiner term
    wire signed [39:0] s1r = accr + (cr[i]*xr[i] - ci[i]*xi[i]);
    wire signed [39:0] s1i = acci + (cr[i]*xi[i] + ci[i]*xr[i]);
    // stage 2: digital combiner term
    wire signed [39:0] s2r = accr + (br[j]*zr[j] - bi[j]*zi[j]);
    wire signed [39:0] s2i = acci + (br[j]*zi[j] + bi[j]*zr[j]);

    always @(posedge clk) begin
        if (rst) begin
            st <= IDLE; cnt <= 0; accr <= 0; acci <= 0; done <= 0;
        end else case (st)
            IDLE: if (start) begin
                st <= S1; cnt <= 0; accr <= 0; acci <= 0; done <= 0;
            end
            S1: begin
                if (i[2:0] == 3'd7) begin           // end of one 8-element sub-array
                    zr[i[5:3]] <= s1r[34:15];
                    zi[i[5:3]] <= s1i[34:15];
                    accr <= 0; acci <= 0;
                end else begin
                    accr <= s1r; acci <= s1i;
                end
                if (cnt == 7'd63) begin st <= S2; cnt <= 0; end
                else cnt <= cnt + 1'b1;
            end
            S2: begin
                accr <= s2r; acci <= s2i;
                if (cnt == 7'd7) st <= FIN;
                else cnt <= cnt + 1'b1;
            end
            FIN: begin
                y_re <= accr[38:15]; y_im <= acci[38:15];
                done <= 1'b1; st <= IDLE;
            end
        endcase
    end
endmodule