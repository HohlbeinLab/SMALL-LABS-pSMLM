function astig_phasor_calibrationRoutine_SL(calibstack,zpositions,savelocation)
[astig_phasor_curveXY,astig_phasor_curveYX] = astig_phasor_calibrationRoutine(calibstack,zpositions);
save(savelocation,'astig_phasor_curveXY','astig_phasor_curveYX')
end