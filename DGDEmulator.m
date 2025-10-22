classdef DGDEmulator
    %CDEMULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static = true)
        function oSigOut = DGDEmulate(oSigIn, DGD, simPara)
            oSigOut = oSigIn;
            oSigOut.compEnvelope = 1/2*TransmissionLink.DGDEmulator.TimeShift(oSigIn.compEnvelope, DGD/2, simPara) +...
                1/2*TransmissionLink.DGDEmulator.TimeShift(oSigIn.compEnvelope, -DGD/2, simPara);            
        end
        function shiftedSampSeq = TimeShift(sampSeq, timeshift, simPara)            
            N = size(sampSeq,1);            
            Omega = 2*pi/N/simPara.samplingPeriod*[0:N/2-1 -N/2:-1].';
            if size(sampSeq,2) == 2
                Omega = [Omega Omega];
            end
            timeShifter = exp(1i*Omega*(timeshift*1e-12));
            sampSeq = ifft(sampSeq);
            sampSeq = timeShifter.*sampSeq;
            shiftedSampSeq = fft(sampSeq);
        end
    end
end

