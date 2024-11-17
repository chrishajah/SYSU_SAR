**T1 Pulse Compression Simulation**
- Requirements
    - 1．LFM信号分析：
        - （1）仿真LFM信号；
        - （2）观察不同过采样率下的DFT结果；
        - （3）观察不同TBP的LFM信号的频谱。
    - 2．脉冲压缩仿真：针对“基带LFM信号”
        - （1）实现无误差的脉冲压缩，计算指标（IRW、PSLR、ISLR）；
        - （2）观察频域加窗的影响，计算指标（IRW、PSLR、ISLR）
    - 3.LFM回波仿真：
        - （1）输入参数:目标参数：RCS等于1,分别位于15Km，16Km。LFM信号参数：中心频率5.0GHz，脉冲宽度30us，带宽30MHz。接收机参数：发射信号后0.1ms启动接收机开始接收信号，持续时间0.04ms，按1.2倍过采样因子采样。
        - （2）输出:完成脉冲压缩。
- Accomplishment
    - 原理：https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/SAR1.pdf
    - 结果：
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/1_1.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/1_2.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/1_3.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/2_1.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/2_2.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/2_3.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/2.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/3_1.png)
    ![img](https://github.com/chrishajah/SYSU_SAR/blob/main/T1_Pulse_Compression/3_2.png)
