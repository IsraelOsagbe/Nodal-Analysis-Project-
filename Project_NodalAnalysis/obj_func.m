function MAPE = obj_func(f1, f2)
    MAPE = (abs((f2 - f1)/f2))*0.5;
end