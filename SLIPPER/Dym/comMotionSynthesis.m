function motionC =  comMotionSynthesis(motionFrame,motionBody,massFrame,massBody)
    motionC = massBody/(massBody+massFrame)*motionBody  +...
              massFrame/(massBody+massFrame)*motionFrame;
end