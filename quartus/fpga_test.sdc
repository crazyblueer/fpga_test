create_clock -name clk -period 20.000 [get_ports clk]
derive_clock_uncertainty
set_false_path -from [get_ports {x_* c_* b_* start rst}]
set_false_path -to   [get_ports {y_* done}]