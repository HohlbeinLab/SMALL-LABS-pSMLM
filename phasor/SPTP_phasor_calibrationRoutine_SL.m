function SPTP_phasor_calibrationRoutine_SL(calibstack,zpositions,savelocation)
[SPTP_phasor_curveX,SPTP_phasor_curveY] = SPTP_phasor_calibrationRoutine(calibstack,zpositions);
save(savelocation,'SPTP_phasor_curveX','SPTP_phasor_curveY')
end