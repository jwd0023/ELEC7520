function vout = R2R(v_reference, ideal_resistance, n_bits, tolerance) 

    % Creates an array of all possible binary inputs for given number of bits.
    binary_inputs = dec2bin([0:(2^(n_bits)-1)])-'0';

    % Create an array of worst-case tolerance values (output resistor (Rout) and 
    % MSB 2R resistor are minimum values, all others are maximums).
    tol_array = ones((n_bits*2)+1,1);
    tol_array(end-1:end) = -tol_array(end-1:end);
    tolerance = tol_array*tolerance;

    % Get vector of ideal resistances (to be modified by the tolerance).
    ideal_resistance = ideal_resistance*ones((n_bits*2)+1,1);

    % Create some 2R values.
    ideal_resistance(1:2) = ideal_resistance(1:2)*2;
    ideal_resistance(4:2:end) = ideal_resistance(4:2:end)*2;

    % Matrix of resistances for ladder network for maximum error.
    resistances = ideal_resistance .* (1+tolerance); 

    % Find base current.
    base_current = v_reference / _find_ladder_req(n_bits, resistances);

    % Find output current.
    output_currents = _find_output_currents(n_bits, v_reference, base_current, resistances);
    output_currents = binary_inputs * output_currents';

    % Compute output voltage at inverting op-amp.
    vout = -output_currents * resistances(end);

end


function req = _find_ladder_req(n_bits, r)

    % Begin with the parallel combination of first two.
    req = 0;
    
    % Build up req progressively through the ladder.
    for i=1:2:(2*n_bits)-1
        req = parallel_resistance([req+r(i), r(i+1)]);
    end
    
end


function resistance = _find_resistance_at_bit(bit, r, n_bits)

    % Begin with the parallel combination of first two.
    denominator = 0;
    
    % Special case if bit == 0.
    if bit ~= 0
        denominator = parallel_resistance([r(1) r(2)]);
    end
    
    % Build up progressively through the ladder.
    for i=3:2:(2*bit)-3
        denominator = parallel_resistance([denominator+r(i), r(i+1)]);
    end
    
    % Always add the last two.
    denominator += r(2*(bit+1)) + r(2*bit + 1);
    
    % Solve final output resistance.
    resistance = r(2*(bit+1))/denominator;
    
end


function output_currents = _find_output_currents(n_bits, v_reference, base_current, r)
    
    % Preallocate output currents.
    output_currents = zeros(1, n_bits);
    
    % Progressively changes based on each current solved.
    total_current = 0;
    
    for i=n_bits-1:-1:0

        resistance = _find_resistance_at_bit(i, r, n_bits);
        output_currents(i+1) = (base_current-total_current) * resistance;
        total_current += output_currents(i+1);
    
    end

    % Flip to match correct bit order.
    output_currents = fliplr(output_currents);
    
end


function req = parallel_resistance(r)
    req = 1.0 ./ sum(1.0 ./ r);
end
