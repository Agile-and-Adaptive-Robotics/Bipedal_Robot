function C = BipedalCost(HumanTorque, RobotTorque, GTorque, GDiameter, diameter, divisions)
%This function is the cost function for the bipedal robot optimization
%calculations. Contributions to the Cost Value:
%   - absolute value differece between the human generated torque and robot
%   - distance between the muscle attacment point and where the bone is
%   - diameter of the festo muscle
%   - length of muscle?


C = 0;

for i = 1:divisions
    for ii = 1:divisions
        C = C + abs(HumanTorque1(i, ii) - RobotTorque1(i, ii));
    end
end

C = C + GDiameter*diameter;