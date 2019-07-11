%________________________________________________________________________%
%1a
%________________________________________________________________________%
%Plots output voltage for an R-2R network and analyzes to see if the
%desired number of bits is possible given the resistor tolerance
Vref = 1; %reference voltage (this has no effect on the final answer, merely scales everything)
Ri = 1; %ideal resistor value (this has no effect on the final answer, merely scales everything)
tol = 0.01; %tolerance
n = 7; %number of bits

vout = R2R(Vref,Ri,n,tol);

voutideal = R2R(Vref,Ri,n,0);

%plot output%
close
figure(2)
hold on
scatter([0:size(vout,1)-1],vout,'filled')
scatter([0:size(vout,1)-1],voutideal)
hold off
title(sprintf("DAC output for n = %d bits and +/- %0.2f%% tolerance\n",[n tol*100]))
xlabel("decimal equivalent of binary number")
ylabel("output voltage")
diff(vout);
if max(diff(vout))>= 0
    fprintf("Monotinicity violated for %d bits and +/- %0.2f%% tolerance\n",[n tol*100])
else
    fprintf("Monotinicity holds for %d bits and +/- %0.2f%% tolerance\n",[n tol*100])
end

function vout = R2R(Vref,Ri,n,tol) 
%computes output voltage of an R-2R network when given the following:
% (1) reference voltage
% (2) design resistor value
% (3) number of bits
% (4) resistor tolerance

%________________________________________________________________________%
%THIS CODE IS ~ALMOST~ ENTIRELY AUTOMATED... NEED TO AUTOMATE SOLUTIONS FOR
%CURRENTS (in) FOR IT TO BE PERFECT, BUT THIS IS CLOSE ENOUGH FOR NOW.. WORKS UP
%TO n=7
%________________________________________________________________________%
b = dec2bin([0:(2^(n)-1)])-'0'; %creates an array of all possible binary inputs for given number of bits

%create an array of worst-case tolerance values (output resistor (Rout) and MSB 2R resistor are minimum values, all others are maximums)
tol_array = ones(n*2+1,1);
tol_array(length(tol_array)-1:length(tol_array)) = -tol_array(length(tol_array)-1:length(tol_array)); %get vector of tolerances
tol = tol_array*tol;


Ri = Ri*ones(n*2+1,1); %get vector of ideal resistances (to be modified by the tolerance)
%2R values are the first resistor and all even numbered array entries
%following
Ri(1:2) = Ri(1:2)*2; %create some 2R values
Ri(4:2:end) = Ri(4:2:end)*2; %create some more 2R values

R = Ri.*(1+tol); %matrix of resistances for ladder network for maximum error


%need to automate this part... solve for all currents in ladder network

%7 bits
if size(b,2) == 7;
    Req = par(par(par(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9),R(10))+R(11),R(12))+R(13),R(14));
    I = Vref/Req;
    itot = 0;
    
    i6 = (I-itot)*R(14)/(par(par(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9),R(10))+R(11),R(12))+R(13)+R(14));
    itot = itot + i6;
    i5 = (I-itot)*R(12)/(par(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9),R(10))+R(11)+R(12));
    itot = itot + i5;
    i4 = (I-itot)*R(10)/(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9)+R(10));
    itot = itot + i4;
    i3 = (I-itot)*R(8)/(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7)+R(8));
    itot = itot + i3;
    i2 = (I-itot)*R(6)/(par(par(R(1),R(2))+R(3),R(4))+R(5)+R(6));
    itot = itot + i2;
    i1 = (I-itot)*R(4)/(par(R(1),R(2))+R(3)+R(4));
    itot = itot + i1;
    i0 = (I-itot)*R(2)/(R(1)+R(2));
    iout = b*[i6;i5;i4;i3;i2;i1;i0];


%6 bits
elseif size(b,2) == 6;
    Req = par(par(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9),R(10))+R(11),R(12));
    I = Vref/Req;
    itot = 0;
    
    i5 = (I-itot)*R(12)/(par(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9),R(10))+R(11)+R(12));
    itot = itot + i5;
    i4 = (I-itot)*R(10)/(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9)+R(10));
    itot = itot + i4;
    i3 = (I-itot)*R(8)/(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7)+R(8));
    itot = itot + i3;
    i2 = (I-itot)*R(6)/(par(par(R(1),R(2))+R(3),R(4))+R(5)+R(6));
    itot = itot + i2;
    i1 = (I-itot)*R(4)/(par(R(1),R(2))+R(3)+R(4));
    itot = itot + i1;
    i0 = (I-itot)*R(2)/(R(1)+R(2));
    iout = b*[i5;i4;i3;i2;i1;i0];


%5 bits
elseif size(b,2) == 5;
    Req = par(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9),R(10));
    I = Vref/Req;
    itot = 0;
    
    i4 = I*R(10)/(par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8))+R(9)+R(10));
    itot = itot + i4;
    i3 = (I-itot)*R(8)/(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7)+R(8));
    itot = itot + i3;
    i2 = (I-itot)*R(6)/(par(par(R(1),R(2))+R(3),R(4))+R(5)+R(6));
    itot = itot + i2;
    i1 = (I-itot)*R(4)/(par(R(1),R(2))+R(3)+R(4));
    itot = itot + i1;
    i0 = (I-itot)*R(2)/(R(1)+R(2));
    iout = b*[i4;i3;i2;i1;i0];


%4 bits
elseif size(b,2) == 4;
    Req = par(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7),R(8));
    I = Vref/Req;
    itot = 0;
        
    i3 = I*R(8)/(par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6))+R(7)+R(8));
    itot = itot + i3;
    i2 = (I-itot)*R(6)/(par(par(R(1),R(2))+R(3),R(4))+R(5)+R(6));
    itot = itot + i2;
    i1 = (I-itot)*R(4)/(par(R(1),R(2))+R(3)+R(4));
    itot = itot + i1;
    i0 = (I-itot)*R(2)/(R(1)+R(2));
    iout = b*[i3;i2;i1;i0];

%3 bits
elseif size(b,2) == 3;
    Req = par(par(par(R(1),R(2))+R(3),R(4))+R(5),R(6));
    I = Vref/Req;
    itot = 0;
    i2 = I*R(6)/(par(par(R(1),R(2))+R(3),R(4))+R(5)+R(6));
    itot = itot + i2;
    i1 = (I-itot)*R(4)/(par(R(1),R(2))+R(3)+R(4));
    itot = itot + i1;
    i0 = (I-itot)*R(2)/(R(1)+R(2));
    iout = b*[i2;i1;i0];

%2 bits
elseif size(b,2) == 2;
    Req = par(par(R(1),R(2))+R(3),R(4));
    I = Vref/Req;
    itot =0;
    i1 = (I-itot)*R(4)/(par(R(1),R(2))+R(3)+R(4));
    itot = itot + i1;
    i0 = (I-itot)*R(2)/(R(1)+R(2));
    iout = b*[i1;i0];

%1 bit
elseif size(b,2) == 1;
    Req = par(R(1),R(2));
    I = Vref/Req;
    i0 = (I)*R(2)/(R(1)+R(2));
    iout = b*[i0];
end

%compute output voltage at inverting op-amp
vout = -iout*R(length(R));

end

function Req = par(R1,R2) %function for combining resistors in parallel
Req = R1*R2/(R1+R2); %if you don't know this you belong in Lowder Hall
end