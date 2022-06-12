/***************************************************
Channel Coding Course Work: conolutional codes
This program template has given the message generator, BPSK modulation, AWGN channel model and BPSK demodulation,
you should first determine the encoder structure, then define the message and codeword length, generate the state table, write the convolutional encoder and decoder.

If you have any question, please contact me via e-mail: wuchy28@mail2.sysu.edu.cn
***************************************************/

#define  _CRT_SECURE_NO_WARNINGS
//#define  DEBUG
//#define  SHOW_message
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<float.h>

#define state_num 2 //the number of the state of encoder structure
#define pow_state_num 4
#define message_length 1024 //the length of message
#define codeword_length 2048 //the length of codeword
#define triple_mes_len 3096	//the length of trubo code
float code_rate = (float)message_length / (float)codeword_length;

// channel coefficient
#define pi 3.1415926
double N0, sgm;

int state_table[pow_state_num][message_length];//state table, the size should be defined yourself
int state_table_systematic[pow_state_num][message_length];//state table, the size should be defined yourself

int message[message_length], codeword[codeword_length];//message and codeword by (7,5)8 encoder
int codeword_systematic[codeword_length];//codeword in systematic way by (1,5/7)8 enocder
int codeword_turbo[triple_mes_len + state_num];//codeword by turbo enocder
int codeword_turbo_puncture[codeword_length + state_num];

int re_codeword[codeword_length], re_codeword_systematic[codeword_length];//the received codeword

int de_message[message_length];//the hard decoding message
int de_message_soft[message_length];//the soft decoding message
int de_message_BCJR[message_length];//the BCJR decoding message
int de_message_systematic[message_length];//the hard decoding message from systematic codeword
int de_message_soft_systematic[message_length];//the soft decoding message from systematic codeword
int de_message_BCJR_systematic[message_length];//the BCJR decoding message from systematic codeword
int de_message_turbo[message_length];// the turbo decoding message

double tx_symbol[codeword_length][2], tx_symbol_systematic[codeword_length][2], tx_symbol_turbo[triple_mes_len + state_num][2], tx_symbol_turbo_puncture[codeword_length + state_num][2];//the transmitted symbols
double rx_symbol[codeword_length][2], rx_symbol_systematic[codeword_length][2], rx_symbol_turbo[triple_mes_len + state_num][2], rx_symbol_turbo_puncture[codeword_length + state_num][2];//the received symbols

double edge_prob[message_length][pow_state_num];//probability of edge in BCJR decoder
double state_former_prob[message_length][pow_state_num];//state probability compute from start
double state_later_prob[message_length][pow_state_num];//state probability compute from end
double posterior_prob[message_length];//the probability of p(x=0|y)

double edge_prob_systematic[message_length][pow_state_num];//probability of edge in BCJR decoder from systematic codeword
double state_former_prob_systematic[message_length][pow_state_num];//state probability compute from start from systematic codeword
double state_later_prob_systematic[message_length][pow_state_num];//state probability compute from end from systematic codeword
double posterior_prob_systematic[message_length];//the probability of p(x=0|y) from systematic codeword

double edge_prob_turbo_p1[message_length][pow_state_num];//probability of edge in BCJR decoder from turbo codeword part 1
double edge_prob_turbo_p2[message_length][pow_state_num];//part 2

void statetable();
void encoder();
void modulation();
void channel();
void demodulation();
void decoder();
void soft_decoder();
int hamming_distance(int *a, int *b, int size);
double Euclidean_distance(double* a, double* b, int size=2);

void BCJR_decoder(double sigma);
void compute_Edge_prob(double sigma);
void compute_State_Former_prob();
void compute_State_Later_prob();
void compute_Posterior_prob();

void encoder_systematic();
void modulation_systematic();
void channel_systematic();
void demodulation_systematic();
void decoder_systematic();
void soft_decoder_systematic();

void BCJR_decoder_systematic(double sigma);
void compute_Edge_prob_systematic(double sigma);
void compute_State_Former_prob_systematic();
void compute_State_Later_prob_systematic();
void compute_Posterior_prob_systematic();

// the main for course work
void main2()
{
	int i;
	float SNR, start, finish;
	long int bit_error, bit_error_soft, bit_error_BCJR, seq, seq_num;
	double BER, BER_soft, BER_BCJR;
	double progress;

	//generate state table
	statetable();

	//random seed
	srand((int)time(0));

	//input the SNR and frame number
	/*printf("\nEnter start SNR: ");
	scanf("%f", &start);
	printf("\nEnter finish SNR: ");
	scanf("%f", &finish);
	printf("\nPlease input the number of message: ");
	scanf("%d", &seq_num);*/
	start = 0;
	finish = 10;
	seq_num = 10000;

	for (SNR = start; SNR <= finish; SNR++)
	{
		//channel noise
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2);

		bit_error = 0;
		bit_error_soft = 0;
		bit_error_BCJR = 0;

		for (seq = 1; seq <= seq_num; seq++)
		{
			//generate binary message randomly
			/****************
			Pay attention that message is appended by 0 whose number is equal to the state of encoder structure.
			****************/

			
			for (i = 0; i < message_length - state_num; i++)
			{
				message[i] = rand() % 2;
			}
			for (i = message_length - state_num; i < message_length; i++)
			{
				message[i] = 0;
			}
			//message[0] = 0; message[1] = 0; message[2] = 0; message[3] = 0; message[4] = 0;
			//message[5] = 0; message[6] = 1; message[7] = 1; message[8] = 0; message[9] = 0;


			//convolutional encoder
			encoder();

			//BPSK modulation
			modulation();

			//AWGN channel
			channel();

			//BPSK demodulation, it's needed in hard-decision Viterbi decoder
			demodulation();

			//convolutional decoder
			decoder();

			//calculate the number of bit error
			for (i = 0; i < message_length; i++)
			{
				if (message[i] != de_message[i])
					bit_error++;
			}

			progress = (double)(seq * 100) / (double)seq_num;

			//calculate the BER
			BER = (double)bit_error / (double)(message_length * seq);

			//****************** for soft decoder ********************
			soft_decoder();

			//calculate the number of bit error
			for (i = 0; i < message_length; i++)
			{
				if (message[i] != de_message_soft[i])
					bit_error_soft++;
			}

			//calculate the BER
			BER_soft = (double)bit_error_soft / (double)(message_length * seq);

			//****************** for BCJR decoder ********************
			BCJR_decoder(sgm);

			//calculate the number of bit error
			for (i = 0; i < message_length; i++)
			{
				if (message[i] != de_message_BCJR[i])
					bit_error_BCJR++;
			}

			//calculate the BER
			BER_BCJR = (double)bit_error_BCJR / (double)(message_length * seq);

			//print the intermediate result
			//printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\r", progress, SNR, bit_error, BER);
			//printf("Progress=%2.1f, SNR=%2.1f, BER=%E, BER_soft=%E, BER_BCJR=%E\r", progress, SNR, BER, BER_soft, BER_BCJR);
		}

		//calculate the BER
		BER = (double)bit_error / (double)(message_length * seq_num);
		BER_soft = (double)bit_error_soft / (double)(message_length * seq_num);
		BER_BCJR = (double)bit_error_BCJR / (double)(message_length * seq_num);

		//print the final result
		//printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n", progress, SNR, bit_error, BER);
		printf("SNR = % 2.1f, BER = % E, BER_soft = % E, BER_BCJR = % E\n", SNR, BER, BER_soft, BER_BCJR);
	}
	system("pause");
}
void statetable()
{
	//allocate the memory fro state_table, a matrix of (pow_state_num * message_length)
	//nothing to do here
}

void encoder()
{
	int current_state = 0;		//initialize state
	//convolution encoder, the input is message[] and the output is codeword[]
	for (int i = 0; i < message_length; ++i) {
		//encode with state transition graph
		if (current_state == 0) {
			if (message[i] == 0) {
				codeword[2 * i] = 0;
				codeword[2 * i + 1] = 0;		//the encode result
				current_state = 0;				//update the next state
			}
			else if (message[i] == 1) {
				codeword[2 * i] = 1;
				codeword[2 * i + 1] = 1;
				current_state = 2;
			}
		}
		else if (current_state == 1) {
			if (message[i] == 0) {
				codeword[2 * i] = 1;
				codeword[2 * i + 1] = 1;
				current_state = 0;
			}
			else if (message[i] == 1) {
				codeword[2 * i] = 0;
				codeword[2 * i + 1] = 0;
				current_state = 2;
			}
		}
		else if (current_state == 2) {
			if (message[i] == 0) {
				codeword[2 * i] = 1;
				codeword[2 * i + 1] = 0;
				current_state = 1;
			}
			else if (message[i] == 1) {
				codeword[2 * i] = 0;
				codeword[2 * i + 1] = 1;
				current_state = 3;
			}
		}
		else if (current_state == 3) {
			if (message[i] == 0) {
				codeword[2 * i] = 0;
				codeword[2 * i + 1] = 1;
				current_state = 1;
			}
			else if (message[i] == 1) {
				codeword[2 * i] = 1;
				codeword[2 * i + 1] = 0;
				current_state = 3;
			}
		}
		
	}
}

void modulation()
{
	//BPSK modulation
	int i;

	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i < codeword_length; i++)
	{
		tx_symbol[i][0] = -1 * (2 * codeword[i] - 1);
		tx_symbol[i][1] = 0;
	}
}
void channel()
{
	//AWGN channel
	int i, j;
	double u, r, g;

	for (i = 0; i < codeword_length; i++)
	{
		for (j = 0; j < 2; j++)
		{
			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));

			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			g = (float)r * cos(2 * pi * u);

			rx_symbol[i][j] = tx_symbol[i][j] + g;
		}
	}
}
void demodulation()
{
	int i;
	double d1, d2;
	for (i = 0; i < codeword_length; i++)
	{
		d1 = (rx_symbol[i][0] - 1) * (rx_symbol[i][0] - 1) + rx_symbol[i][1] * rx_symbol[i][1];
		d2 = (rx_symbol[i][0] + 1) * (rx_symbol[i][0] + 1) + rx_symbol[i][1] * rx_symbol[i][1];
		if (d1 < d2)
			re_codeword[i] = 0;
		else
			re_codeword[i] = 1;
	}
}

void decoder()
{
	int state_hamming_distance[4] = { 100000,10000,10000,10000 };		//set the hamming distance for each state large
	state_hamming_distance[0] = 0;		//state initialize to 0, so the hamming distance of state 0 is 0
	int a[2] = { 0,0 };
	int b[2] = { 0,1 };
	int c[2] = { 1,0 };
	int d[2] = { 1,1 };					//for computing hamming distance between re_codeword and encode output for specific state

	for (int i = 0; i < message_length; ++i) {
		int mes_part[2];
		mes_part[0] = re_codeword[2 * i];
		mes_part[1] = re_codeword[2 * i + 1];		//fetch a pair of deocde input
		int tmp1, tmp2, s0, s1, s2, s3;

		//update state_hamming_distance and state_table
		tmp1 = state_hamming_distance[0] + hamming_distance(a, mes_part, state_num);
		tmp2 = state_hamming_distance[1] + hamming_distance(d, mes_part, state_num);		
			//for state 0 in next state, compute two possible routine and select the routine with minimum hamming distance
		state_table[0][i] = tmp1 < tmp2 ? 0 : 1;		//store the former state number, that means routine
		s0 = tmp1 < tmp2 ? tmp1 : tmp2;					//the hamming value of state 0 to be update

		tmp1 = state_hamming_distance[2] + hamming_distance(c, mes_part, state_num);
		tmp2 = state_hamming_distance[3] + hamming_distance(b, mes_part, state_num);
		state_table[1][i] =  tmp1 < tmp2 ? 2 : 3;
		s1 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_hamming_distance[0] + hamming_distance(d, mes_part, state_num);
		tmp2 = state_hamming_distance[1] + hamming_distance(a, mes_part, state_num);
		state_table[2][i] = tmp1 < tmp2 ? 0 : 1;
		s2 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_hamming_distance[2] + hamming_distance(b, mes_part, state_num);
		tmp2 = state_hamming_distance[3] + hamming_distance(c, mes_part, state_num);
		state_table[3][i] = tmp1 < tmp2 ? 2 : 3;
		s3 = tmp1 < tmp2 ? tmp1 : tmp2;

		state_hamming_distance[0] = s0;
		state_hamming_distance[1] = s1;
		state_hamming_distance[2] = s2;
		state_hamming_distance[3] = s3;					//update hamming value for 4 state
		/*
		printf("state_hamming_distance:%d,%d,%d,%d \n", state_hamming_distance[0], state_hamming_distance[1], \
			state_hamming_distance[2], state_hamming_distance[3]);
			*/
	}


	//decode in reverse way
	int reverse_state = 0;		//final state is 0
	for (int i = message_length - 1; i >= 0; i--) {
		switch (reverse_state) {
		case 0:de_message[i] = 0; break;
		case 1:de_message[i] = 0; break;		//given the later state, message is decisive
		case 2:de_message[i] = 1; break;
		case 3:de_message[i] = 1; break;
		}
		reverse_state = state_table[reverse_state][i];		//update the former state, from end to start
	}
}
int hamming_distance(int* a, int* b, int size) {		//compute hamming distance of two vectors
	int result = 0;
	for (int i = 0; i < size; ++i) {
		result += (a[i] == b[i] ? 0 : 1);				//element not equal, hamming distance + 1
	}
	return result;
}

void soft_decoder()
{

	double state_Euclidean_distance[4] = { 100000,10000,10000,10000 };		//set the Euclidean distance for each state large
	state_Euclidean_distance[0] = 0;		//state initialize to 0, so the Euclidean distance of state 0 is 0
	double a[2] = { 1,0 };
	double b[2] = { -1,0 };				//for computing Euclidean distance between rx_symbol and modulation output for specific state

	for (int i = 0; i < message_length; ++i) {
		double tmp1, tmp2, s0, s1, s2, s3;

		//update state_Euclidean_distance and state_table
		tmp1 = state_Euclidean_distance[0] + Euclidean_distance(a, rx_symbol[2 * i]) + Euclidean_distance(a, rx_symbol[2 * i + 1]);
		tmp2 = state_Euclidean_distance[1] + Euclidean_distance(b, rx_symbol[2 * i]) + Euclidean_distance(b, rx_symbol[2 * i + 1]);
		//for state 0 in next state, compute two possible routine and select the routine with minimum Euclidean distance
		state_table[0][i] = tmp1 < tmp2 ? 0 : 1;		//store the former state number, that means routine
		s0 = tmp1 < tmp2 ? tmp1 : tmp2;					//the Euclidean value of state 0 to be update

		tmp1 = state_Euclidean_distance[2] + Euclidean_distance(b, rx_symbol[2 * i]) + Euclidean_distance(a, rx_symbol[2 * i + 1]);
		tmp2 = state_Euclidean_distance[3] + Euclidean_distance(a, rx_symbol[2 * i]) + Euclidean_distance(b, rx_symbol[2 * i + 1]);
		state_table[1][i] = tmp1 < tmp2 ? 2 : 3;
		s1 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_Euclidean_distance[0] + Euclidean_distance(b, rx_symbol[2 * i]) + Euclidean_distance(b, rx_symbol[2 * i + 1]);
		tmp2 = state_Euclidean_distance[1] + Euclidean_distance(a, rx_symbol[2 * i]) + Euclidean_distance(a, rx_symbol[2 * i + 1]);
		state_table[2][i] = tmp1 < tmp2 ? 0 : 1;
		s2 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_Euclidean_distance[2] + Euclidean_distance(a, rx_symbol[2 * i]) + Euclidean_distance(b, rx_symbol[2 * i + 1]);//hamming_distance(b, mes_part, state_num);
		tmp2 = state_Euclidean_distance[3] + Euclidean_distance(b, rx_symbol[2 * i]) + Euclidean_distance(a, rx_symbol[2 * i + 1]);//hamming_distance(c, mes_part, state_num);
		state_table[3][i] = tmp1 < tmp2 ? 2 : 3;
		s3 = tmp1 < tmp2 ? tmp1 : tmp2;

		state_Euclidean_distance[0] = s0;
		state_Euclidean_distance[1] = s1;
		state_Euclidean_distance[2] = s2;
		state_Euclidean_distance[3] = s3;					//update Euclidean value for 4 state
	}


	//decode in reverse way
	int reverse_state = 0;		//final state is 0
	for (int i = message_length - 1; i >= 0; i--) {
		switch (reverse_state) {
		case 0:de_message_soft[i] = 0; break;
		case 1:de_message_soft[i] = 0; break;		//given the later state, message is decisive
		case 2:de_message_soft[i] = 1; break;
		case 3:de_message_soft[i] = 1; break;
		}
		reverse_state = state_table[reverse_state][i];		//update the former state, from end to start
	}
}
double Euclidean_distance(double* a, double* b, int size)
{
	double result = 0;
	for (int i = 0; i < size; ++i) {
		result += (a[i]-b[i])*(a[i]-b[i]);				//Euclidean distance between vector a and b
	}
	return result;
}

void BCJR_decoder(double sigma)
{
	compute_Edge_prob(sigma);		// need to know noise power, which is a problem
	compute_State_Former_prob();	// compute state probability from start to end
	compute_State_Later_prob();		// compute state probability from end to start
	compute_Posterior_prob();		// compute p(x=0|y)

	for (int i = 0; i < message_length; ++i) {	// the decode logic
		if (posterior_prob[i] >= 0.5)
			de_message_BCJR[i] = 0;				// if p(x=0|y)>0.5 then decode the symbol as 0. else as 1
		else
			de_message_BCJR[i] = 1;
	}
}
void compute_Edge_prob(double sigma)
{
	double p0[state_num];
	double p1[state_num];		// state_num indicate how many receieve symbols the edge probability related to
	for (int i = 0; i < message_length; ++i) {
		for (int j = 0; j < state_num; ++j) {		// compute and store channel observasion for an edge probability
			p0[j] = exp(-(rx_symbol[i * 2 + j][0] - 1) * (rx_symbol[i * 2 + j][0] - 1) / 2 / sigma / sigma);	//probability of receive symbol to be 0 according to channel observation
			p1[j] = exp(-(rx_symbol[i * 2 + j][0] + 1) * (rx_symbol[i * 2 + j][0] + 1) / 2 / sigma / sigma);	//				...					  1					...
			p0[j] = p0[j] / (p0[j] + p1[j]);		//set sum of prob. is 1
			p1[j] = 1 - p0[j];
		}

		edge_prob[i][0] = p0[0] * p0[1];
		edge_prob[i][1] = p0[0] * p1[1];
		edge_prob[i][2] = p1[0] * p0[1];
		edge_prob[i][3] = p1[0] * p1[1];		//the edge probability, with no information about x, hence p(x=0)=p(x=1)=0.5
	}
}
void compute_State_Former_prob()
{
	double p0, p1, p2, p3, p_sum;
	state_former_prob[0][0] = 1;
	state_former_prob[0][1] = 0;
	state_former_prob[0][2] = 0;
	state_former_prob[0][3] = 0;				//the very start of state probability, with 00 state has probability 1 and other have 0
	for (int i = 1; i < message_length; ++i) {
		p0 = state_former_prob[i - 1][0] * edge_prob[i - 1][0] + state_former_prob[i - 1][1] * edge_prob[i - 1][3];
		p1 = state_former_prob[i - 1][2] * edge_prob[i - 1][2] + state_former_prob[i - 1][3] * edge_prob[i - 1][1];
		p2 = state_former_prob[i - 1][0] * edge_prob[i - 1][3] + state_former_prob[i - 1][1] * edge_prob[i - 1][0];
		p3 = state_former_prob[i - 1][2] * edge_prob[i - 1][1] + state_former_prob[i - 1][3] * edge_prob[i - 1][2];		//compute probability of next state according to trellis

		p_sum = p0 + p1 + p2 + p3;
		state_former_prob[i][0] = p0 / p_sum;
		state_former_prob[i][1] = p1 / p_sum;
		state_former_prob[i][2] = p2 / p_sum;
		state_former_prob[i][3] = p3 / p_sum;		//normalize the probability and store into matrix
	}
}
void compute_State_Later_prob()
{
	double p0, p1, p2, p3, p_sum;
	state_later_prob[message_length - 1][0] = 1;
	state_later_prob[message_length - 1][1] = 0;
	state_later_prob[message_length - 1][2] = 0;
	state_later_prob[message_length - 1][3] = 0;	//the very end of state probability, with 00 state has probability 1 and other have 0
	for (int i = message_length - 2; i >= 0; --i) {
		p0 = state_later_prob[i + 1][0] * edge_prob[i + 1][0] + state_later_prob[i + 1][2] * edge_prob[i + 1][3];
		p1 = state_later_prob[i + 1][0] * edge_prob[i + 1][3] + state_later_prob[i + 1][2] * edge_prob[i + 1][0];
		p2 = state_later_prob[i + 1][1] * edge_prob[i + 1][2] + state_later_prob[i + 1][3] * edge_prob[i + 1][1];
		p3 = state_later_prob[i + 1][1] * edge_prob[i + 1][1] + state_later_prob[i + 1][3] * edge_prob[i + 1][2];		//compute probability of former state according to trellis

		p_sum = p0 + p1 + p2 + p3;
		state_later_prob[i][0] = p0 / p_sum;
		state_later_prob[i][1] = p1 / p_sum;
		state_later_prob[i][2] = p2 / p_sum;
		state_later_prob[i][3] = p3 / p_sum;		//normalize the probability and store into matrix
	}
}
void compute_Posterior_prob()
{
	double p0, p1;
	for (int i = 0; i < message_length; ++i) {
		p0 = state_former_prob[i][0] * edge_prob[i][0] * state_later_prob[i][0] + \
			state_former_prob[i][1] * edge_prob[i][3] * state_later_prob[i][0] + \
			state_former_prob[i][2] * edge_prob[i][2] * state_later_prob[i][1] + \
			state_former_prob[i][3] * edge_prob[i][1] * state_later_prob[i][1];			//every [start status -- edge -- end status] that cause by an input 0

		p1 = state_former_prob[i][0] * edge_prob[i][3] * state_later_prob[i][2] + \
			state_former_prob[i][1] * edge_prob[i][0] * state_later_prob[i][2] + \
			state_former_prob[i][2] * edge_prob[i][1] * state_later_prob[i][3] + \
			state_former_prob[i][3] * edge_prob[i][2] * state_later_prob[i][3];			//every [start status -- edge -- end status] that cause by an input 1

		posterior_prob[i] = p0 / (p0 + p1);				//the probability of decode output 0
	}
}

/* --------the coursework code is done above------- */

/// <summary>
/// next is the systematic implementation for fun
/// </summary>

void encoder_systematic()
{
	int current_state = 0;		//initialize state
	//convolution encoder, the input is message[] and the output is codeword_systematic[]
	//note that the last 2 bits in message are not always 00, so discard the last 2 bits
	for (int i = 0; i < message_length-2; ++i) {
		//encode with state transition graph
		if (current_state == 0) {
			if (message[i] == 0) {
				codeword_systematic[2 * i] = 0;
				codeword_systematic[2 * i + 1] = 0;		//the encode result
				current_state = 0;				//update the next state
			}
			else if (message[i] == 1) {
				codeword_systematic[2 * i] = 1;
				codeword_systematic[2 * i + 1] = 1;
				current_state = 2;
			}
		}
		else if (current_state == 1) {
			if (message[i] == 0) {
				codeword_systematic[2 * i] = 0;
				codeword_systematic[2 * i + 1] = 0;
				current_state = 2;
			}
			else if (message[i] == 1) {
				codeword_systematic[2 * i] = 1;
				codeword_systematic[2 * i + 1] = 1;
				current_state = 0;
			}
		}
		else if (current_state == 2) {
			if (message[i] == 0) {
				codeword_systematic[2 * i] = 0;
				codeword_systematic[2 * i + 1] = 1;
				current_state = 3;
			}
			else if (message[i] == 1) {
				codeword_systematic[2 * i] = 1;
				codeword_systematic[2 * i + 1] = 0;
				current_state = 1;
			}
		}
		else if (current_state == 3) {
			if (message[i] == 0) {
				codeword_systematic[2 * i] = 0;
				codeword_systematic[2 * i + 1] = 1;
				current_state = 1;
			}
			else if (message[i] == 1) {
				codeword_systematic[2 * i] = 1;
				codeword_systematic[2 * i + 1] = 0;
				current_state = 3;
			}
		}
	}
	//add the last 2 bits in message to let state be 0
	if (current_state == 0) {
		codeword_systematic[2 * (message_length - 2)] = 0;
		codeword_systematic[2 * (message_length - 2) + 1] = 0;
		codeword_systematic[2 * (message_length - 1)] = 0;
		codeword_systematic[2 * (message_length - 1) + 1] = 0;
	}
	else if (current_state == 1) {
		codeword_systematic[2 * (message_length - 2)] = 1;
		codeword_systematic[2 * (message_length - 2) + 1] = 1;
		codeword_systematic[2 * (message_length - 1)] = 0;
		codeword_systematic[2 * (message_length - 1) + 1] = 0;
	}
	else if (current_state == 2) {
		codeword_systematic[2 * (message_length - 2)] = 1;
		codeword_systematic[2 * (message_length - 2) + 1] = 0;
		codeword_systematic[2 * (message_length - 1)] = 1;
		codeword_systematic[2 * (message_length - 1) + 1] = 1;
	}
	else {
		codeword_systematic[2 * (message_length - 2)] = 0;
		codeword_systematic[2 * (message_length - 2) + 1] = 1;
		codeword_systematic[2 * (message_length - 1)] = 1;
		codeword_systematic[2 * (message_length - 1) + 1] = 1;
	}
}

void modulation_systematic()
{
	//BPSK modulation
	int i;

	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i < codeword_length; i++)
	{
		tx_symbol_systematic[i][0] = (double)-1 * (2 * codeword_systematic[i] - 1);
		tx_symbol_systematic[i][1] = 0;
	}
}
void channel_systematic()
{
	//AWGN channel
	int i, j;
	double u, r, g;

	for (i = 0; i < codeword_length; i++)
	{
		for (j = 0; j < 2; j++)
		{
			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));

			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			g = (float)r * cos(2 * pi * u);

			rx_symbol_systematic[i][j] = tx_symbol_systematic[i][j] + g;
		}
	}
}
void demodulation_systematic()
{
	int i;
	double d1, d2;
	for (i = 0; i < codeword_length; i++)
	{
		d1 = (rx_symbol_systematic[i][0] - 1) * (rx_symbol_systematic[i][0] - 1) + rx_symbol_systematic[i][1] * rx_symbol_systematic[i][1];
		d2 = (rx_symbol_systematic[i][0] + 1) * (rx_symbol_systematic[i][0] + 1) + rx_symbol_systematic[i][1] * rx_symbol_systematic[i][1];
		if (d1 < d2)
			re_codeword_systematic[i] = 0;
		else
			re_codeword_systematic[i] = 1;
	}
}

void decoder_systematic() {
	int state_hamming_distance[4] = { 100000,10000,10000,10000 };		//set the hamming distance for each state large
	state_hamming_distance[0] = 0;		//state initialize to 0, so the hamming distance of state 0 is 0
	int a[2] = { 0,0 };
	int b[2] = { 0,1 };
	int c[2] = { 1,0 };
	int d[2] = { 1,1 };					//for computing hamming distance between re_codeword and encode output for specific state

	for (int i = 0; i < message_length; ++i) {
		int mes_part[2];
		mes_part[0] = re_codeword_systematic[2 * i];
		mes_part[1] = re_codeword_systematic[2 * i + 1];		//fetch a pair of deocde input
		int tmp1, tmp2, s0, s1, s2, s3;

		//update state_hamming_distance and state_table
		tmp1 = state_hamming_distance[0] + hamming_distance(a, mes_part, state_num);
		tmp2 = state_hamming_distance[1] + hamming_distance(d, mes_part, state_num);
		//for state 0 in next state, compute two possible routine and select the routine with minimum hamming distance
		state_table_systematic[0][i] = tmp1 < tmp2 ? 0 : 1;		//store the former state number, that means routine
		s0 = tmp1 < tmp2 ? tmp1 : tmp2;					//the hamming value of state 0 to be update

		tmp1 = state_hamming_distance[2] + hamming_distance(c, mes_part, state_num);
		tmp2 = state_hamming_distance[3] + hamming_distance(b, mes_part, state_num);
		state_table_systematic[1][i] = tmp1 < tmp2 ? 2 : 3;
		s1 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_hamming_distance[0] + hamming_distance(d, mes_part, state_num);
		tmp2 = state_hamming_distance[1] + hamming_distance(a, mes_part, state_num);
		state_table_systematic[2][i] = tmp1 < tmp2 ? 0 : 1;
		s2 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_hamming_distance[2] + hamming_distance(b, mes_part, state_num);
		tmp2 = state_hamming_distance[3] + hamming_distance(c, mes_part, state_num);
		state_table_systematic[3][i] = tmp1 < tmp2 ? 2 : 3;
		s3 = tmp1 < tmp2 ? tmp1 : tmp2;

		state_hamming_distance[0] = s0;
		state_hamming_distance[1] = s1;
		state_hamming_distance[2] = s2;
		state_hamming_distance[3] = s3;					//update hamming value for 4 state
		/*
		printf("state_hamming_distance:%d,%d,%d,%d \n", state_hamming_distance[0], state_hamming_distance[1], \
			state_hamming_distance[2], state_hamming_distance[3]);
			*/
	}
	//decode in reverse way
	int reverse_state = 0;		//final state is 0
	int former_reverse_state;
	for (int i = message_length - 1; i >= 0; i--) {
		former_reverse_state = state_table_systematic[reverse_state][i];		//update the former state, from end to start
		switch (reverse_state) {
		case 0:de_message_systematic[i] = (former_reverse_state == 0 ? 0 : 1); break;	//the decode output depend on the later state and former state
		case 1:de_message_systematic[i] = (former_reverse_state == 2 ? 1 : 0); break;
		case 2:de_message_systematic[i] = (former_reverse_state == 0 ? 1 : 0); break;
		case 3:de_message_systematic[i] = (former_reverse_state == 2 ? 0 : 1); break;
		}
		reverse_state = former_reverse_state;		//update the former state, from end to start
	}
}
void soft_decoder_systematic()
{
	double state_Euclidean_distance[4] = { 100000,10000,10000,10000 };		//set the Euclidean distance for each state large
	state_Euclidean_distance[0] = 0;		//state initialize to 0, so the Euclidean distance of state 0 is 0
	double a[2] = { 1,0 };
	double b[2] = { -1,0 };				//for computing Euclidean distance between rx_symbol and modulation output for specific state

	for (int i = 0; i < message_length; ++i) {
		double tmp1, tmp2, s0, s1, s2, s3;

		//update state_Euclidean_distance and state_table
		tmp1 = state_Euclidean_distance[0] + Euclidean_distance(a, rx_symbol_systematic[2 * i]) + Euclidean_distance(a, rx_symbol_systematic[2 * i + 1]);
		tmp2 = state_Euclidean_distance[1] + Euclidean_distance(b, rx_symbol_systematic[2 * i]) + Euclidean_distance(b, rx_symbol_systematic[2 * i + 1]);
		//for state 0 in next state, compute two possible routine and select the routine with minimum Euclidean distance
		state_table_systematic[0][i] = tmp1 < tmp2 ? 0 : 1;		//store the former state number, that means routine
		s0 = tmp1 < tmp2 ? tmp1 : tmp2;					//the Euclidean value of state 0 to be update

		tmp1 = state_Euclidean_distance[2] + Euclidean_distance(b, rx_symbol_systematic[2 * i]) + Euclidean_distance(a, rx_symbol_systematic[2 * i + 1]);
		tmp2 = state_Euclidean_distance[3] + Euclidean_distance(a, rx_symbol_systematic[2 * i]) + Euclidean_distance(b, rx_symbol_systematic[2 * i + 1]);
		state_table_systematic[1][i] = tmp1 < tmp2 ? 2 : 3;
		s1 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_Euclidean_distance[0] + Euclidean_distance(b, rx_symbol_systematic[2 * i]) + Euclidean_distance(b, rx_symbol_systematic[2 * i + 1]);
		tmp2 = state_Euclidean_distance[1] + Euclidean_distance(a, rx_symbol_systematic[2 * i]) + Euclidean_distance(a, rx_symbol_systematic[2 * i + 1]);
		state_table_systematic[2][i] = tmp1 < tmp2 ? 0 : 1;
		s2 = tmp1 < tmp2 ? tmp1 : tmp2;

		tmp1 = state_Euclidean_distance[2] + Euclidean_distance(a, rx_symbol_systematic[2 * i]) + Euclidean_distance(b, rx_symbol_systematic[2 * i + 1]);//hamming_distance(b, mes_part, state_num);
		tmp2 = state_Euclidean_distance[3] + Euclidean_distance(b, rx_symbol_systematic[2 * i]) + Euclidean_distance(a, rx_symbol_systematic[2 * i + 1]);//hamming_distance(c, mes_part, state_num);
		state_table_systematic[3][i] = tmp1 < tmp2 ? 2 : 3;
		s3 = tmp1 < tmp2 ? tmp1 : tmp2;

		state_Euclidean_distance[0] = s0;
		state_Euclidean_distance[1] = s1;
		state_Euclidean_distance[2] = s2;
		state_Euclidean_distance[3] = s3;					//update Euclidean value for 4 state
	}


	//decode in reverse way
	int reverse_state = 0;		//final state is 0
	int former_reverse_state;
	for (int i = message_length - 1; i >= 0; i--) {
		former_reverse_state = state_table_systematic[reverse_state][i];		//update the former state, from end to start
		switch (reverse_state) {
		case 0:de_message_soft_systematic[i] = (former_reverse_state == 0 ? 0 : 1); break;	//the decode output depend on the later state and former state
		case 1:de_message_soft_systematic[i] = (former_reverse_state == 2 ? 1 : 0); break;
		case 2:de_message_soft_systematic[i] = (former_reverse_state == 0 ? 1 : 0); break;
		case 3:de_message_soft_systematic[i] = (former_reverse_state == 2 ? 0 : 1); break;
		}
		reverse_state = former_reverse_state;		//update the former state, from end to start
	}
}

void BCJR_decoder_systematic(double sigma)
{
	compute_Edge_prob_systematic(sigma);		// need to know noise power, which is a problem
	compute_State_Former_prob_systematic();		// compute state probability from start to end
	compute_State_Later_prob_systematic();		// compute state probability from end to start
	compute_Posterior_prob_systematic();		// compute p(x=0|y)

	for (int i = 0; i < message_length; ++i) {	// the decode logic
		if (posterior_prob_systematic[i] >= 0.5)
			de_message_BCJR_systematic[i] = 0;				// if p(x=0|y)>0.5 then decode the symbol as 0. else as 1
		else
			de_message_BCJR_systematic[i] = 1;
	}
}
void compute_Edge_prob_systematic(double sigma)
{
	double p0[state_num];
	double p1[state_num];		// state_num indicate how many receieve symbols the edge probability related to
	for (int i = 0; i < message_length; ++i) {
		for (int j = 0; j < state_num; ++j) {		// compute and store channel observasion for an edge probability
			p0[j] = exp(-(rx_symbol_systematic[i * 2 + j][0] - 1) * (rx_symbol_systematic[i * 2 + j][0] - 1) / 2 / sigma / sigma);	//probability of receive symbol to be 0 according to channel observation
			p1[j] = exp(-(rx_symbol_systematic[i * 2 + j][0] + 1) * (rx_symbol_systematic[i * 2 + j][0] + 1) / 2 / sigma / sigma);	//				...					  1					...
			p0[j] = p0[j] / (p0[j] + p1[j]);		//set sum of prob. is 1
			p1[j] = 1 - p0[j];
		}

		edge_prob_systematic[i][0] = p0[0] * p0[1];
		edge_prob_systematic[i][1] = p0[0] * p1[1];
		edge_prob_systematic[i][2] = p1[0] * p0[1];
		edge_prob_systematic[i][3] = p1[0] * p1[1];		//the edge probability, with no information about x, hence p(x=0)=p(x=1)=0.5
	}
}
void compute_State_Former_prob_systematic()
{
	double p0, p1, p2, p3, p_sum;
	state_former_prob_systematic[0][0] = 1;
	state_former_prob_systematic[0][1] = 0;
	state_former_prob_systematic[0][2] = 0;
	state_former_prob_systematic[0][3] = 0;				//the very start of state probability, with 00 state has probability 1 and other have 0
	for (int i = 1; i < message_length; ++i) {
		p0 = state_former_prob_systematic[i - 1][0] * edge_prob_systematic[i - 1][0] + state_former_prob_systematic[i - 1][1] * edge_prob_systematic[i - 1][3];
		p1 = state_former_prob_systematic[i - 1][2] * edge_prob_systematic[i - 1][2] + state_former_prob_systematic[i - 1][3] * edge_prob_systematic[i - 1][1];
		p2 = state_former_prob_systematic[i - 1][0] * edge_prob_systematic[i - 1][3] + state_former_prob_systematic[i - 1][1] * edge_prob_systematic[i - 1][0];
		p3 = state_former_prob_systematic[i - 1][2] * edge_prob_systematic[i - 1][1] + state_former_prob_systematic[i - 1][3] * edge_prob_systematic[i - 1][2];		//compute probability of next state according to trellis

		p_sum = p0 + p1 + p2 + p3;
		state_former_prob_systematic[i][0] = p0 / p_sum;
		state_former_prob_systematic[i][1] = p1 / p_sum;
		state_former_prob_systematic[i][2] = p2 / p_sum;
		state_former_prob_systematic[i][3] = p3 / p_sum;		//normalize the probability and store into matrix
	}	
}
void compute_State_Later_prob_systematic()
{
	double p0, p1, p2, p3, p_sum;
	state_later_prob_systematic[message_length - 1][0] = 1;
	state_later_prob_systematic[message_length - 1][1] = 0;
	state_later_prob_systematic[message_length - 1][2] = 0;
	state_later_prob_systematic[message_length - 1][3] = 0;	//the very end of state probability, with 00 state has probability 1 and other have 0
	for (int i = message_length - 2; i >= 0; --i) {
		p0 = state_later_prob_systematic[i + 1][0] * edge_prob_systematic[i + 1][0] + state_later_prob_systematic[i + 1][2] * edge_prob_systematic[i + 1][3];
		p1 = state_later_prob_systematic[i + 1][0] * edge_prob_systematic[i + 1][3] + state_later_prob_systematic[i + 1][2] * edge_prob_systematic[i + 1][0];
		p2 = state_later_prob_systematic[i + 1][1] * edge_prob_systematic[i + 1][2] + state_later_prob_systematic[i + 1][3] * edge_prob_systematic[i + 1][1];
		p3 = state_later_prob_systematic[i + 1][1] * edge_prob_systematic[i + 1][1] + state_later_prob_systematic[i + 1][3] * edge_prob_systematic[i + 1][2];		//compute probability of former state according to trellis

		p_sum = p0 + p1 + p2 + p3;
		state_later_prob_systematic[i][0] = p0 / p_sum;
		state_later_prob_systematic[i][1] = p1 / p_sum;
		state_later_prob_systematic[i][2] = p2 / p_sum;
		state_later_prob_systematic[i][3] = p3 / p_sum;		//normalize the probability and store into matrix
	}
}
void compute_Posterior_prob_systematic()
{
	double p0, p1;
	for (int i = 0; i < message_length; ++i) {
		p0 = state_former_prob_systematic[i][0] * edge_prob_systematic[i][0] * state_later_prob_systematic[i][0] + \
			state_former_prob_systematic[i][1] * edge_prob_systematic[i][0] * state_later_prob_systematic[i][2] + \
			state_former_prob_systematic[i][2] * edge_prob_systematic[i][1] * state_later_prob_systematic[i][3] + \
			state_former_prob_systematic[i][3] * edge_prob_systematic[i][1] * state_later_prob_systematic[i][1];			//every [start status -- edge -- end status] that cause by an input 0

		p1 = state_former_prob_systematic[i][0] * edge_prob_systematic[i][3] * state_later_prob_systematic[i][2] + \
			state_former_prob_systematic[i][1] * edge_prob_systematic[i][3] * state_later_prob_systematic[i][0] + \
			state_former_prob_systematic[i][2] * edge_prob_systematic[i][2] * state_later_prob_systematic[i][1] + \
			state_former_prob_systematic[i][3] * edge_prob_systematic[i][2] * state_later_prob_systematic[i][3];			//every [start status -- edge -- end status] that cause by an input 1

		posterior_prob_systematic[i] = p0 / (p0 + p1);				//the probability of decode output 0
	}
}

/// <summary>
/// next is the turbo code implementation for fun
/// </summary>

void interleave(const int* message, int* interleaved_message);		// interleave message, where interleaved_message should be allocated
void de_interleave(const int* interleaved_message, int* message);		// de-interleave message
void interleave(const double* message, double* interleaved_message);		// interleave message, where interleaved_message should be allocated
void de_interleave(const double* interleaved_message, double* message);	// de-interleave message
void encoder_systematic_siso(const int* message, int* encoded_message);	// systematic conv.encoder without systematic output
void encoder_turbo();				// turbo encoder with code rate 1/3
void puncture();					// puncture encoded turbo code, setting coderate as 1/2

void modulation_turbo();
void channel_turbo();
void modulation_turbo_puncture();
void channel_turbo_puncture();
void restore_puncture();

void decoder_turbo(const double sigma, const int iteration=10);								//turbo decode function
void compute_Edge_prob_turbo(const double sigma);										//compute edge prob with out Px, which can be done out of the iteration
void BCJR_decoder_turbo_p1(const double* Px/* prior prob. P(x=0)*/, double* Px_y/* posterior prob. P(x=0|y)*/, double* Pex/* exterior prob. Pe(x=0)*/);
void BCJR_decoder_turbo_p2(const double* Px/* prior prob. P(x=0)*/, double* Px_y/* posterior prob. P(x=0|y)*/, double* Pex/* exterior prob. Pe(x=0)*/);
void compute_State_Former_prob_trubo_p1(const double* Px/* prior prob. */);
void compute_State_Former_prob_trubo_p2(const double* Px/* prior prob. */);
void compute_State_Later_prob_trubo_p1(const double* Px/* prior prob. */);
void compute_State_Later_prob_trubo_p2(const double* Px/* prior prob. */);
void compute_Posterior_prob_trubo_p1(const double* Px/* prior prob. */);
void compute_Posterior_prob_trubo_p2(const double* Px/* prior prob. */);			//counterpart to systematic but with prior prob. p(x)
double prevent_0(double input);

int random_interleave_coefficient[message_length - 2] =
{ 408, 869, 942, 168, 500, 335, 107, 760, 401, 564, 242, 950, 296, 995, 212, 764, 162, 198, 120, 880, 639, 92, 274, 245, 835, 551, 361, 885, 585, 32, 939, 724, 973, 809, 320, 906, 479, 321, 1006, 579, 700, 608, 222, 169, 65, 260, 630, 66, 478, 466, 114, 917, 907, 794, 64, 559, 912, 915, 153, 90, 205, 295, 317, 392, 77, 124, 790, 944, 558, 404, 499, 5, 951, 81, 727, 220, 480, 193, 132, 359, 619, 683, 145, 228, 922, 291, 409, 1004, 956, 286, 712, 675, 269, 792, 831, 221, 758, 482, 701, 464, 781, 190, 483, 1009, 855, 101, 3, 412, 14, 537, 780, 142, 121, 612, 225, 504, 836, 518, 259, 15, 50, 444, 195, 416, 959, 104, 717, 626, 488, 757, 1007, 272, 767, 916, 397, 173, 271, 474, 946, 67, 713, 870, 845, 775, 495, 388, 958, 419, 834, 22, 661, 862, 277, 438, 659, 723, 38, 645, 961, 770, 829, 599, 571, 473, 514, 1015, 768, 215, 591, 603, 20, 878, 657, 117, 56, 247, 233, 89, 976, 888, 71, 27, 207, 131, 873, 230, 600, 106, 625, 678, 583, 343, 606, 227, 549, 652, 761, 246, 133, 797, 270, 415, 688, 351, 847, 379, 341, 879, 743, 128, 919, 570, 577, 12, 735, 428, 642, 432, 934, 979, 468, 986, 992, 360, 347, 772, 33, 704, 530, 816, 545, 346, 467, 560, 447, 298, 348, 620, 237, 884, 1000, 714, 968, 574, 111, 436, 384, 337, 325, 858, 202, 99, 926, 364, 807, 708, 550, 867, 676, 505, 344, 87, 414, 214, 139, 967, 533, 947, 618, 46, 801, 183, 413, 234, 471, 842, 36, 256, 838, 670, 266, 850, 352, 410, 982, 720, 1012, 1003, 354, 658, 394, 7, 693, 771, 521, 899, 178, 283, 422, 261, 990, 170, 93, 433, 1020, 425, 696, 581, 244, 345, 112, 522, 75, 427, 350, 788, 17, 595, 923, 324, 769, 726, 886, 980, 161, 371, 681, 1013, 985, 526, 918, 539, 928, 598, 13, 426, 251, 974, 866, 588, 323, 48, 110, 957, 59, 510, 154, 52, 602, 892, 647, 857, 865, 975, 634, 330, 509, 689, 116, 800, 812, 332, 991, 821, 507, 703, 860, 44, 26, 403, 531, 186, 200, 282, 452, 555, 1017, 725, 853, 502, 1011, 863, 840, 340, 442, 83, 54, 279, 833, 152, 548, 138, 1021, 825, 804, 1016, 903, 406, 86, 441, 322, 782, 481, 126, 313, 305, 163, 34, 827, 63, 287, 894, 664, 257, 232, 55, 308, 799, 590, 617, 963, 424, 458, 216, 924, 293, 765, 381, 143, 294, 70, 830, 492, 752, 783, 370, 336, 806, 503, 311, 334, 57, 594, 188, 776, 498, 451, 728, 516, 303, 1002, 872, 43, 949, 16, 965, 301, 597, 938, 534, 338, 935, 715, 380, 402, 589, 310, 766, 543, 893, 276, 472, 557, 365, 856, 908, 289, 476, 368, 469, 751, 952, 484, 339, 326, 805, 243, 118, 964, 23, 777, 837, 791, 541, 122, 580, 937, 576, 418, 74, 265, 4, 672, 508, 72, 278, 29, 280, 832, 1018, 73, 546, 366, 268, 970, 137, 687, 465, 709, 569, 211, 876, 561, 113, 275, 231, 1001, 669, 515, 349, 948, 496, 30, 88, 258, 686, 648, 568, 206, 849, 249, 930, 774, 523, 181, 638, 622, 184, 798, 252, 802, 219, 864, 297, 213, 966, 706, 996, 999, 927, 641, 637, 994, 881, 160, 778, 998, 874, 674, 1014, 632, 646, 196, 604, 129, 889, 435, 182, 817, 695, 229, 385, 898, 921, 905, 592, 68, 140, 95, 931, 819, 40, 981, 11, 902, 264, 175, 21, 563, 871, 654, 936, 374, 281, 746, 315, 431, 440, 997, 954, 1010, 10, 859, 628, 31, 115, 826, 983, 528, 624, 382, 762, 123, 164, 37, 58, 151, 682, 623, 0, 61, 201, 644, 319, 532, 665, 192, 203, 738, 993, 731, 627, 396, 455, 690, 813, 750, 828, 610, 239, 165, 248, 984, 158, 566, 9, 943, 290, 653, 578, 517, 611, 240, 971, 736, 895, 722, 393, 929, 877, 565, 130, 582, 953, 1005, 640, 756, 556, 742, 463, 744, 42, 519, 485, 204, 718, 475, 136, 353, 573, 633, 697, 405, 820, 24, 988, 443, 395, 386, 721, 538, 839, 629, 667, 39, 972, 144, 684, 607, 155, 562, 786, 814, 933, 430, 377, 312, 49, 784, 940, 489, 199, 304, 176, 85, 2, 314, 470, 400, 962, 486, 491, 96, 376, 734, 740, 650, 171, 739, 177, 179, 745, 146, 94, 875, 694, 818, 19, 527, 823, 614, 824, 60, 852, 300, 6, 166, 896, 41, 900, 180, 754, 197, 909, 223, 375, 891, 329, 25, 803, 53, 796, 529, 977, 76, 631, 267, 378, 407, 284, 808, 391, 224, 209, 662, 82, 680, 733, 613, 789, 914, 255, 449, 535, 932, 217, 253, 318, 331, 453, 461, 666, 945, 411, 141, 730, 811, 45, 651, 513, 747, 460, 333, 732, 553, 226, 494, 848, 901, 911, 157, 439, 1, 655, 910, 854, 759, 98, 446, 134, 238, 729, 102, 355, 208, 299, 660, 635, 506, 292, 367, 779, 194, 429, 913, 621, 810, 846, 103, 544, 685, 80, 149, 501, 390, 119, 679, 462, 707, 373, 448, 601, 596, 91, 512, 51, 955, 716, 147, 941, 649, 263, 69, 369, 711, 785, 357, 35, 159, 445, 719, 841, 47, 572, 108, 882, 897, 587, 191, 586, 524, 493, 434, 887, 8, 851, 540, 189, 307, 218, 389, 421, 174, 372, 542, 254, 236, 815, 960, 477, 520, 605, 698, 691, 309, 844, 748, 288, 925, 363, 773, 185, 868, 135, 459, 795, 62, 793, 399, 616, 387, 609, 456, 362, 100, 710, 398, 741, 969, 904, 692, 890, 127, 84, 536, 987, 861, 358, 668, 250, 273, 187, 702, 677, 18, 417, 172, 167, 525, 547, 420, 643, 78, 327, 306, 584, 79, 450, 497, 673, 763, 671, 328, 316, 749, 28, 148, 487, 656, 787, 125, 342, 755, 920, 1019, 663, 636, 737, 262, 210, 356, 978, 150, 615, 423, 575, 552, 457, 753, 1008, 285, 454, 235, 105, 97, 593, 156, 109, 822, 699, 511, 241, 705, 567, 843, 989, 302, 383, 883, 490, 437, 554 };
//{ 8, 49, 23, 67, 77, 61, 76, 46, 7, 4, 19, 72, 45, 30, 94, 71, 81, 35, 87, 64, 25, 56, 41, 44, 82, 21, 84, 40, 33, 79, 26, 32, 78, 52, 39, 0, 18, 93, 88, 27, 16, 15, 34, 58, 24, 38, 3, 13, 57, 70, 12, 36, 17, 90, 53, 1, 69, 55, 50, 83, 62, 91, 68, 75, 43, 80, 97, 96, 2, 47, 92, 29, 89, 60, 22, 37, 31, 20, 54, 66, 28, 85, 42, 6, 48, 14, 10, 65, 63, 95, 74, 5, 86, 59, 51, 11, 9, 73 };

void interleave(const int* message, int* interleaved_message)
{
	for (int i = 0; i < message_length - state_num; ++i) {
		interleaved_message[i] = message[random_interleave_coefficient[i]];
	}
	for (int i = 0; i < state_num; ++i) {
		interleaved_message[message_length - 2 + i] = message[message_length - 2 + i];
	}
}
void de_interleave(const int* interleaved_message, int* message)
{
	for (int i = 0; i < message_length - state_num; ++i) {
		message[random_interleave_coefficient[i]] = interleaved_message[i];
	}
	for (int i = 0; i < state_num; ++i) {
		message[message_length - 2 + i] = interleaved_message[message_length - 2 + i];
	}
}
void interleave(const double* message, double* interleaved_message)
{
	for (int i = 0; i < message_length - state_num; ++i) {
		interleaved_message[i] = message[random_interleave_coefficient[i]];
	}
	for (int i = 0; i < state_num; ++i) {
		interleaved_message[message_length - 2 + i] = message[message_length - 2 + i];
	}
}
void de_interleave(const double* interleaved_message, double* message)
{
	for (int i = 0; i < message_length - state_num; ++i) {
		message[random_interleave_coefficient[i]] = interleaved_message[i];
	}
	for (int i = 0; i < state_num; ++i) {
		message[message_length - 2 + i] = interleaved_message[message_length - 2 + i];
	}
}
void encoder_systematic_siso(const int* message, int* encoded_message)
{
	int current_state = 0;		//initialize state
	//convolution encoder, the input is message[] and the output is encoded_message[]
	//note that the last 2 bits in message are not always 00, so discard the last 2 bits
	for (int i = 0; i < message_length - 2; ++i) {
		//encode with state transition graph
		if (current_state == 0) {
			if (message[i] == 0) {
				encoded_message[i] = 0;		//the encode result
				current_state = 0;				//update the next state
			}
			else if (message[i] == 1) {
				encoded_message[i] = 1;
				current_state = 2;
			}
		}
		else if (current_state == 1) {
			if (message[i] == 0) {
				encoded_message[i] = 0;
				current_state = 2;
			}
			else if (message[i] == 1) {
				encoded_message[i] = 1;
				current_state = 0;
			}
		}
		else if (current_state == 2) {
			if (message[i] == 0) {
				encoded_message[i] = 1;
				current_state = 3;
			}
			else if (message[i] == 1) {
				encoded_message[i] = 0;
				current_state = 1;
			}
		}
		else if (current_state == 3) {
			if (message[i] == 0) {
				encoded_message[i] = 1;
				current_state = 1;
			}
			else if (message[i] == 1) {
				encoded_message[i] = 0;
				current_state = 3;
			}
		}
	}
	//add the last 2 bits in message to let state be 0
	if (current_state == 0) {
		encoded_message[(message_length - 2)] = 0;
		encoded_message[(message_length - 1)] = 0;
		// ------------- the later two bits are the correspoding systematic last two end bits ----------- 
		encoded_message[(message_length)] = 0;
		encoded_message[(message_length + 1)] = 0;
	}
	else if (current_state == 1) {
		encoded_message[(message_length - 2)] = 1;
		encoded_message[(message_length - 1)] = 0;
		encoded_message[(message_length)] = 1;
		encoded_message[(message_length + 1)] = 0;
	}
	else if (current_state == 2) {
		encoded_message[(message_length - 2)] = 0;
		encoded_message[(message_length - 1)] = 1;
		encoded_message[(message_length)] = 1;
		encoded_message[(message_length + 1)] = 1;
	}
	else {
		encoded_message[(message_length - 2)] = 1;
		encoded_message[(message_length - 1)] = 1;
		encoded_message[(message_length)] = 0;
		encoded_message[(message_length + 1)] = 1;
	}
}
void encoder_turbo()
{
	int tmp1[message_length + state_num], tmp2[message_length + state_num];
	encoder_systematic_siso(message, tmp1);		//first encoder output
	int interleaved_message[message_length];
	interleave(message, interleaved_message);	//interleave the message
	encoder_systematic_siso(interleaved_message, tmp2);		//second encoder output

	for (int i = 0; i < message_length; ++i) {
		codeword_turbo[3 * i] = message[i];		//systematic part of codeword
		codeword_turbo[3 * i + 1] = tmp1[i];	//non-systematic part 1 of codeword
		codeword_turbo[3 * i + 2] = tmp2[i];	//non-systematic part 2 of codeword
	}

	//make the codeword_turbo include the last systematic bits
	codeword_turbo[3 * (message_length - 2)] = tmp1[message_length];
	codeword_turbo[3 * (message_length - 1)] = tmp1[message_length + 1];
	codeword_turbo[triple_mes_len] = tmp2[message_length];
	codeword_turbo[triple_mes_len + 1] = tmp2[message_length + 1];
}
void puncture()
{
	// puncture patern : [1 0; 0 1]
	for (int i = 0; i < message_length; ++i) {
		codeword_turbo_puncture[2 * i] = codeword_turbo[3 * i];							//the systematic part
		codeword_turbo_puncture[2 * i + 1] = codeword_turbo[3 * i + 1 + (i % 2 == 1)];		//the non-systematic part
	}
	codeword_turbo_puncture[message_length] = codeword_turbo[triple_mes_len];
	codeword_turbo_puncture[message_length + 1] = codeword_turbo[triple_mes_len + 1];		//the systematic part for part2 of turbo code, add to the end of codeword
}

void modulation_turbo()
{
	//BPSK modulation
	int i;

	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i < triple_mes_len + state_num; i++)
	{
		tx_symbol_turbo[i][0] = -1 * (2 * codeword_turbo[i] - 1);
		tx_symbol_turbo[i][1] = 0;
	}
}
void channel_turbo()
{
	//AWGN channel
	int i, j;
	double u, r, g;

	for (i = 0; i < triple_mes_len + state_num; i++)
	{
		for (j = 0; j < 2; j++)
		{
			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));

			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			g = (float)r * cos(2 * pi * u);

			rx_symbol_turbo[i][j] = tx_symbol_turbo[i][j] + g;
		}
	}
}
void modulation_turbo_puncture()
{
	//BPSK modulation
	int i;

	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i < codeword_length + state_num; i++)
	{
		tx_symbol_turbo_puncture[i][0] = -1 * (2 * codeword_turbo_puncture[i] - 1);
		tx_symbol_turbo_puncture[i][1] = 0;
	}
}
void channel_turbo_puncture()
{
	//AWGN channel
	int i, j;
	double u, r, g;

	for (i = 0; i < codeword_length + state_num; i++)
	{
		for (j = 0; j < 2; j++)
		{
			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));

			u = (float)rand() / (float)RAND_MAX;
			if (u == 1.0)
				u = 0.999999;
			g = (float)r * cos(2 * pi * u);

			rx_symbol_turbo_puncture[i][j] = tx_symbol_turbo_puncture[i][j] + g;
		}
	}
}
void restore_puncture() 
{
	//turn rx_symbol_turbo_puncture into rx_symbol_turbo for decoding
	for (int i = 0; i < message_length; ++i) {
		for (int j = 0; j < 2; ++j) {				//NOTE: this decode vector is shared with turbo decoder vector
			rx_symbol_turbo[3 * i][j] = rx_symbol_turbo_puncture[2 * i][j];		//the systematic part is transmitted and recieved
			rx_symbol_turbo[3 * i + 1][j] = (i % 2 == 0) ? rx_symbol_turbo_puncture[2 * i + 1][j] : 0;	//non-systematic part
			rx_symbol_turbo[3 * i + 2][j] = (i % 2 == 1) ? rx_symbol_turbo_puncture[2 * i + 1][j] : 0;
		}
	}
	//the last two systematic bits for part2
	for (int j = 0; j < 2; ++j) {
		rx_symbol_turbo[triple_mes_len][j] = rx_symbol_turbo_puncture[codeword_length][j];
		rx_symbol_turbo[triple_mes_len + 1][j] = rx_symbol_turbo_puncture[codeword_length + 1][j];
	}
}

void compute_Edge_prob_turbo(const double sigma)
{
	double sigma2;
	if (sigma < 0.1)
		sigma2 = 0.1;
	else
		sigma2 = sigma;
	//first generate the interleaved rx_symbol_turbo systemetic part
	double tmp[message_length];		//store rx_symbol_turbo systemetic part as an vector
	for (int i = 0; i < message_length - 2; ++i) {
		tmp[i] = rx_symbol_turbo[3 * i][0];
	}
	tmp[message_length - 2] = rx_symbol_turbo[triple_mes_len][0];
	tmp[message_length - 1] = rx_symbol_turbo[triple_mes_len + 1][0];
	double rx_symbol_I_turbo_interleaved[message_length];	//the I part for interleaved systematic turbo code
	interleave(tmp, rx_symbol_I_turbo_interleaved);

	double p0[state_num * 2];		// times 2 means compute two parts edge prob in turbo
	double p1[state_num * 2];		// state_num indicate how many receieve symbols the edge probability related to
	for (int i = 0; i < message_length; ++i) {
		// compute and store channel observasion for an edge probability

		//for systematic part 1
		p0[0] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_turbo[3 * i][0] - 1) * (rx_symbol_turbo[3 * i][0] - 1) / 2 / sigma2 / sigma2);	//probability of receive symbol to be 0 according to channel observation
		p1[0] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_turbo[3 * i][0] + 1) * (rx_symbol_turbo[3 * i][0] + 1) / 2 / sigma2 / sigma2);	//				...					  1					...
		p0[0] = prevent_0(p0[0]) / (prevent_0(p0[0]) + prevent_0(p1[0]));		//set sum of prob. is 1
		p1[0] = 1 - p0[0];

		//for non-systematic part 1
		p0[1] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_turbo[3 * i + 1][0] - 1) * (rx_symbol_turbo[3 * i + 1][0] - 1) / 2 / sigma2 / sigma2);	//probability of receive symbol to be 0 according to channel observation
		p1[1] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_turbo[3 * i + 1][0] + 1) * (rx_symbol_turbo[3 * i + 1][0] + 1) / 2 / sigma2 / sigma2);	//				...					  1					...
		p0[1] = prevent_0(p0[1]) / (prevent_0(p0[1]) + prevent_0(p1[1]));		//set sum of prob. is 1
		p1[1] = 1 - p0[1];

		//for systematic part 2
		p0[2] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_I_turbo_interleaved[i] - 1) * (rx_symbol_I_turbo_interleaved[i] - 1) / 2 / sigma2 / sigma2);	//probability of receive symbol to be 0 according to channel observation
		p1[2] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_I_turbo_interleaved[i] + 1) * (rx_symbol_I_turbo_interleaved[i] + 1) / 2 / sigma2 / sigma2);	//				...					  1					...
		p0[2] = prevent_0(p0[2]) / (prevent_0(p0[2]) + prevent_0(p1[2]));		//set sum of prob. is 1
		p1[2] = 1 - p0[2];

		//for non-systematic part 2
		p0[3] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_turbo[3 * i + 2][0] - 1) * (rx_symbol_turbo[3 * i + 2][0] - 1) / 2 / sigma2 / sigma2);	//probability of receive symbol to be 0 according to channel observation
		p1[3] = rx_symbol_turbo[3 * i + 2][0] == 0 ? 0.5 : exp(-(rx_symbol_turbo[3 * i + 2][0] + 1) * (rx_symbol_turbo[3 * i + 2][0] + 1) / 2 / sigma2 / sigma2);	//				...					  1					...
		p0[3] = prevent_0(p0[3]) / (prevent_0(p0[3]) + prevent_0(p1[3]));		//set sum of prob. is 1
		p1[3] = 1 - p0[3];

		//compute edge_prob
		edge_prob_turbo_p1[i][0] = p0[0] * p0[1];		//NOTE: do not share the variable with systematic BCJR decoder
		edge_prob_turbo_p1[i][1] = p0[0] * p1[1];
		edge_prob_turbo_p1[i][2] = p1[0] * p0[1];
		edge_prob_turbo_p1[i][3] = p1[0] * p1[1];		//the edge probability, with no information about x, hence p(x=0)=p(x=1)=0.5

		edge_prob_turbo_p2[i][0] = p0[2] * p0[3];		//NOTE: do not share the variable with systematic BCJR decoder
		edge_prob_turbo_p2[i][1] = p0[2] * p1[3];
		edge_prob_turbo_p2[i][2] = p1[2] * p0[3];
		edge_prob_turbo_p2[i][3] = p1[2] * p1[3];		//the edge probability, with no information about x, hence p(x=0)=p(x=1)=0.5
#ifdef DEBUG
		printf("p2[%d]:%f %f %f %f\n", i, p0[2], p0[3], p1[2], p1[3]);
		printf("p1[%d]:%f %f %f %f; ", i,edge_prob_turbo_p1[i][0], edge_prob_turbo_p1[i][1], edge_prob_turbo_p1[i][2], edge_prob_turbo_p1[i][3]);
		printf("p2[%d]:%f %f %f %f\n", i,edge_prob_turbo_p2[i][0], edge_prob_turbo_p2[i][1], edge_prob_turbo_p2[i][2], edge_prob_turbo_p2[i][3]);
#endif	//DEBUG
	}
}
void BCJR_decoder_turbo_p1(const double* Px/* prior prob. */, double* Px_y/* posterior prob. */, double* Pex/* exterior prob. */)
{
	compute_State_Former_prob_trubo_p1(Px/* prior prob. */);		// compute state probability from start to end
	compute_State_Later_prob_trubo_p1(Px/* prior prob. */);		// compute state probability from end to start
	compute_Posterior_prob_trubo_p1(Px/* prior prob. */);		// compute p(x=0|y)

	for (int i = 0; i < message_length; ++i) {	// the decode logic
		Px_y[i] = posterior_prob_systematic[i];
		Pex[i] = Px_y[i] / prevent_0(Px[i]) / (Px_y[i] / prevent_0(Px[i]) + (1 - Px_y[i]) / prevent_0(1 - Px[i]));// compute and normalize
#ifdef DEBUG
		printf("[p1] Px[%d]:%f, ", i, Px[i]);
		printf("Px_y[%d]:%f, ", i, Px_y[i]);
		printf("new Pex[%d]:%f \n", i, Pex[i]);
#endif // DEBUG		
	}
}
void compute_State_Former_prob_trubo_p1(const double* Px/* prior prob. */)
{
	double p0, p1, p2, p3, p_sum;
	state_former_prob_systematic[0][0] = 1;				//NOTE: share the variable with systematic decoder
	state_former_prob_systematic[0][1] = 0;
	state_former_prob_systematic[0][2] = 0;
	state_former_prob_systematic[0][3] = 0;				//the very start of state probability, with 00 state has probability 1 and other have 0
	for (int i = 1; i < message_length; ++i) {
		p0 = state_former_prob_systematic[i - 1][0] * edge_prob_turbo_p1[i - 1][0] * Px[i - 1] + state_former_prob_systematic[i - 1][1] * edge_prob_turbo_p1[i - 1][3] * (1 - Px[i - 1]);
		p1 = state_former_prob_systematic[i - 1][2] * edge_prob_turbo_p1[i - 1][2] * (1 - Px[i - 1]) + state_former_prob_systematic[i - 1][3] * edge_prob_turbo_p1[i - 1][1] * Px[i - 1];
		p2 = state_former_prob_systematic[i - 1][0] * edge_prob_turbo_p1[i - 1][3] * (1 - Px[i - 1]) + state_former_prob_systematic[i - 1][1] * edge_prob_turbo_p1[i - 1][0] * Px[i - 1];
		p3 = state_former_prob_systematic[i - 1][2] * edge_prob_turbo_p1[i - 1][1] * Px[i - 1] + state_former_prob_systematic[i - 1][3] * edge_prob_turbo_p1[i - 1][2] * (1 - Px[i - 1]);		//compute probability of next state according to trellis

		p_sum = prevent_0(p0 + (p1) + (p2) + (p3));
		state_former_prob_systematic[i][0] = (p0) / p_sum;
		state_former_prob_systematic[i][1] = (p1) / p_sum;
		state_former_prob_systematic[i][2] = (p2) / p_sum;
		state_former_prob_systematic[i][3] = (p3) / p_sum;		//normalize the probability and store into matrix
#ifdef DEBUG
		printf("[p1] p0:%f, p1:%f, p2:%f, p3:%f, p_sum:%f\n", p0, p1, p2, p3, p_sum);
		printf("[p1] state_former_prob_systematic[%d]:%f, %f, %f, %f; Px:%f \n", i, state_former_prob_systematic[i][0], state_former_prob_systematic[i][1], state_former_prob_systematic[i][2], state_former_prob_systematic[i][3], Px[i - 1]);
#endif // DEBUG
	}
}
void compute_State_Later_prob_trubo_p1(const double* Px/* prior prob. */)
{
	double p0, p1, p2, p3, p_sum;
	state_later_prob_systematic[message_length - 1][0] = 1;	//NOTE: share the variable with systematic decoder
	state_later_prob_systematic[message_length - 1][1] = 0;
	state_later_prob_systematic[message_length - 1][2] = 0;
	state_later_prob_systematic[message_length - 1][3] = 0;	//the very end of state probability, with 00 state has probability 1 and other have 0
	for (int i = message_length - 2; i >= 0; --i) {
		p0 = state_later_prob_systematic[i + 1][0] * edge_prob_turbo_p1[i + 1][0] * Px[i + 1] + state_later_prob_systematic[i + 1][2] * edge_prob_turbo_p1[i + 1][3] * (1 - Px[i + 1]);
		p1 = state_later_prob_systematic[i + 1][0] * edge_prob_turbo_p1[i + 1][3] * (1 - Px[i + 1]) + state_later_prob_systematic[i + 1][2] * edge_prob_turbo_p1[i + 1][0] * Px[i + 1];
		p2 = state_later_prob_systematic[i + 1][1] * edge_prob_turbo_p1[i + 1][2] * (1 - Px[i + 1]) + state_later_prob_systematic[i + 1][3] * edge_prob_turbo_p1[i + 1][1] * Px[i + 1];
		p3 = state_later_prob_systematic[i + 1][1] * edge_prob_turbo_p1[i + 1][1] * Px[i + 1] + state_later_prob_systematic[i + 1][3] * edge_prob_turbo_p1[i + 1][2] * (1 - Px[i + 1]);		//compute probability of former state according to trellis

		p_sum = prevent_0(p0 + p1 + p2 + p3);
		state_later_prob_systematic[i][0] = (p0) / p_sum;
		state_later_prob_systematic[i][1] = (p1) / p_sum;
		state_later_prob_systematic[i][2] = (p2) / p_sum;
		state_later_prob_systematic[i][3] = (p3) / p_sum;		//normalize the probability and store into matrix
#ifdef DEBUG
		printf("[p1] p0:%f, p1:%f, p2:%f, p3:%f, p_sum:%f\n", p0, p1, p2, p3, p_sum);
		printf("[p1] state_later_prob_systematic[%d]:%f, %f, %f, %f; Px:%f \n", i, state_later_prob_systematic[i][0], state_later_prob_systematic[i][1], state_later_prob_systematic[i][2], state_later_prob_systematic[i][3], Px[i + 1]);
#endif // DEBUG
	}
}
void compute_Posterior_prob_trubo_p1(const double* Px/* prior prob. */)
{
	double p0, p1;
	for (int i = 0; i < message_length; ++i) {
		p0 = state_former_prob_systematic[i][0] * edge_prob_turbo_p1[i][0] * state_later_prob_systematic[i][0] + \
			state_former_prob_systematic[i][1] * edge_prob_turbo_p1[i][0] * state_later_prob_systematic[i][2] + \
			state_former_prob_systematic[i][2] * edge_prob_turbo_p1[i][1] * state_later_prob_systematic[i][3] + \
			state_former_prob_systematic[i][3] * edge_prob_turbo_p1[i][1] * state_later_prob_systematic[i][1];			//every [start status -- edge -- end status] that cause by an input 0

		p1 = state_former_prob_systematic[i][0] * edge_prob_turbo_p1[i][3] * state_later_prob_systematic[i][2] + \
			state_former_prob_systematic[i][1] * edge_prob_turbo_p1[i][3] * state_later_prob_systematic[i][0] + \
			state_former_prob_systematic[i][2] * edge_prob_turbo_p1[i][2] * state_later_prob_systematic[i][1] + \
			state_former_prob_systematic[i][3] * edge_prob_turbo_p1[i][2] * state_later_prob_systematic[i][3];			//every [start status -- edge -- end status] that cause by an input 1
		
		posterior_prob_systematic[i] = Px[i] * prevent_0(p0) / (Px[i] * prevent_0(p0) + (1 - Px[i]) * prevent_0(p1));				//the probability of decode output 0
#ifdef DEBUG
		printf("p0:%f ", p0);
		printf("p1:%f \n", p1);
#endif // DEBUG
	}
}

void BCJR_decoder_turbo_p2(const double* Px/* prior prob. */, double* Px_y/* posterior prob. */, double* Pex/* exterior prob. */)
{
	compute_State_Former_prob_trubo_p2(Px/* prior prob. */);		// compute state probability from start to end
	compute_State_Later_prob_trubo_p2(Px/* prior prob. */);		// compute state probability from end to start
	compute_Posterior_prob_trubo_p2(Px/* prior prob. */);		// compute p(x=0|y)

	for (int i = 0; i < message_length; ++i) {	// the decode logic
		Px_y[i] = posterior_prob_systematic[i];
		Pex[i] = Px_y[i] / prevent_0(Px[i]) / (Px_y[i] / prevent_0(Px[i]) + (1 - Px_y[i]) / prevent_0(1 - Px[i]));// compute and normalize
#ifdef DEBUG
		printf("[p2] Px[%d]:%f, ", i, Px[i]);
		printf("Px_y[%d]:%f, ", i, Px_y[i]);
		printf("Pex[%d]:%f \n", i, Pex[i]);
#endif // DEBUG
	}
}
void compute_State_Former_prob_trubo_p2(const double* Px/* prior prob. */)
{
	double p0, p1, p2, p3, p_sum;
	state_former_prob_systematic[0][0] = 1;				//NOTE: share the variable with systematic decoder
	state_former_prob_systematic[0][1] = 0;
	state_former_prob_systematic[0][2] = 0;
	state_former_prob_systematic[0][3] = 0;				//the very start of state probability, with 00 state has probability 1 and other have 0
	for (int i = 1; i < message_length; ++i) {
		p0 = state_former_prob_systematic[i - 1][0] * edge_prob_turbo_p2[i - 1][0] * Px[i - 1] + state_former_prob_systematic[i - 1][1] * edge_prob_turbo_p2[i - 1][3] * (1 - Px[i - 1]);
		p1 = state_former_prob_systematic[i - 1][2] * edge_prob_turbo_p2[i - 1][2] * (1 - Px[i - 1]) + state_former_prob_systematic[i - 1][3] * edge_prob_turbo_p2[i - 1][1] * Px[i - 1];
		p2 = state_former_prob_systematic[i - 1][0] * edge_prob_turbo_p2[i - 1][3] * (1 - Px[i - 1]) + state_former_prob_systematic[i - 1][1] * edge_prob_turbo_p2[i - 1][0] * Px[i - 1];
		p3 = state_former_prob_systematic[i - 1][2] * edge_prob_turbo_p2[i - 1][1] * Px[i - 1] + state_former_prob_systematic[i - 1][3] * edge_prob_turbo_p2[i - 1][2] * (1 - Px[i - 1]);		//compute probability of next state according to trellis

		p_sum = prevent_0(p0 + (p1) + (p2)+ (p3));
		state_former_prob_systematic[i][0] = (p0) / p_sum;
		state_former_prob_systematic[i][1] = (p1) / p_sum;
		state_former_prob_systematic[i][2] = (p2) / p_sum;
		state_former_prob_systematic[i][3] = (p3) / p_sum;		//normalize the probability and store into matrix
#ifdef DEBUG
		printf("[p2] state_former_prob_systematic[%d]:%f, %f, %f, %f; Px:%f \n", i, state_former_prob_systematic[i][0], state_former_prob_systematic[i][1], state_former_prob_systematic[i][2], state_former_prob_systematic[i][3], Px[i - 1]);
#endif // DEBUG

	}
}
void compute_State_Later_prob_trubo_p2(const double* Px/* prior prob. */)
{
	double p0, p1, p2, p3, p_sum;
	state_later_prob_systematic[message_length - 1][0] = 1;	//NOTE: share the variable with systematic decoder
	state_later_prob_systematic[message_length - 1][1] = 0;
	state_later_prob_systematic[message_length - 1][2] = 0;
	state_later_prob_systematic[message_length - 1][3] = 0;	//the very end of state probability, with 00 state has probability 1 and other have 0
	for (int i = message_length - 2; i >= 0; --i) {
		p0 = state_later_prob_systematic[i + 1][0] * edge_prob_turbo_p2[i + 1][0] * Px[i + 1] + state_later_prob_systematic[i + 1][2] * edge_prob_turbo_p2[i + 1][3] * (1 - Px[i + 1]);
		p1 = state_later_prob_systematic[i + 1][0] * edge_prob_turbo_p2[i + 1][3] * (1 - Px[i + 1]) + state_later_prob_systematic[i + 1][2] * edge_prob_turbo_p2[i + 1][0] * Px[i + 1];
		p2 = state_later_prob_systematic[i + 1][1] * edge_prob_turbo_p2[i + 1][2] * (1 - Px[i + 1]) + state_later_prob_systematic[i + 1][3] * edge_prob_turbo_p2[i + 1][1] * Px[i + 1];
		p3 = state_later_prob_systematic[i + 1][1] * edge_prob_turbo_p2[i + 1][1] * Px[i + 1] + state_later_prob_systematic[i + 1][3] * edge_prob_turbo_p2[i + 1][2] * (1 - Px[i + 1]);		//compute probability of former state according to trellis

		p_sum = prevent_0(p0 + (p1) + (p2) + (p3));
		state_later_prob_systematic[i][0] = (p0) / p_sum;
		state_later_prob_systematic[i][1] = (p1) / p_sum;
		state_later_prob_systematic[i][2] = (p2) / p_sum;
		state_later_prob_systematic[i][3] = (p3) / p_sum;		//normalize the probability and store into matrix
#ifdef DEBUG
		printf("[p2] p0:%f, p1:%f, p2:%f, p3:%f, p_sum:%f\n", p0, p1, p2, p3, p_sum);
		printf("[p2] state_later_prob_systematic[%d]:%f, %f, %f, %f; Px:%f \n", i, state_later_prob_systematic[i][0], state_later_prob_systematic[i][1], state_later_prob_systematic[i][2], state_later_prob_systematic[i][3], Px[i + 1]);
#endif // DEBUG
	}
}
void compute_Posterior_prob_trubo_p2(const double* Px/* prior prob. */)
{
	double p0, p1;
	for (int i = 0; i < message_length; ++i) {
		p0 = state_former_prob_systematic[i][0] * edge_prob_turbo_p2[i][0] * state_later_prob_systematic[i][0] + \
			state_former_prob_systematic[i][1] * edge_prob_turbo_p2[i][0] * state_later_prob_systematic[i][2] + \
			state_former_prob_systematic[i][2] * edge_prob_turbo_p2[i][1] * state_later_prob_systematic[i][3] + \
			state_former_prob_systematic[i][3] * edge_prob_turbo_p2[i][1] * state_later_prob_systematic[i][1];			//every [start status -- edge -- end status] that cause by an input 0

		p1 = state_former_prob_systematic[i][0] * edge_prob_turbo_p2[i][3] * state_later_prob_systematic[i][2] + \
			state_former_prob_systematic[i][1] * edge_prob_turbo_p2[i][3] * state_later_prob_systematic[i][0] + \
			state_former_prob_systematic[i][2] * edge_prob_turbo_p2[i][2] * state_later_prob_systematic[i][1] + \
			state_former_prob_systematic[i][3] * edge_prob_turbo_p2[i][2] * state_later_prob_systematic[i][3];			//every [start status -- edge -- end status] that cause by an input 1
#ifdef DEBUG
		printf("p0:%f ", p0);
		printf("p1:%f \n", p1);
#endif // DEBUG

		posterior_prob_systematic[i] = Px[i] * prevent_0(p0) / (Px[i] * prevent_0(p0) + (1 - Px[i]) * prevent_0(p1));				//the probability of decode output 0
	}
}
double prevent_0(double input)
{
	if (input == 0) {
		return DBL_MIN;
	}
	else
		return input;
}

void decoder_turbo(const double sigma, const int iteration)
{
	compute_Edge_prob_turbo(sigma);
	double Px[message_length];
	for (int i = 0; i < message_length; ++i) {
		Px[i] = 0.5;
	}
	double Px_y[message_length];
	double Pex[message_length];
	double p1_Pex_last[state_num];
	double p2_Pex_last[state_num];			//the last two bits of posterior probability need to be stored specifically
	for (int j = 0; j < state_num; ++j) {
		p1_Pex_last[j] = 0.5;
		p2_Pex_last[j] = 0.5;				//at first we don't know the systematic prob.
	}
	for (int i = 0; i < iteration; ++i) {
		BCJR_decoder_turbo_p1(Px, Px_y, Pex);
		for (int j = 0; j < state_num; ++j) {
			p1_Pex_last[j] = Pex[message_length - 2 + j];
		}
		interleave(Pex, Px);
		for (int j = 0; j < state_num; ++j) {
			Px[message_length - 2 + j] = p2_Pex_last[j];
		}
		BCJR_decoder_turbo_p2(Px, Px_y, Pex);
		for (int j = 0; j < state_num; ++j) {
			p2_Pex_last[j] = Pex[message_length - 2 + j];
		}
		de_interleave(Pex, Px);
		for (int j = 0; j < state_num; ++j) {
			Px[message_length - 2 + j] = p1_Pex_last[j];
		}
#ifdef DEBUG
		printf("\n!!!!!\n");
		for (int j = 0; j < message_length; ++j) {
			printf("%f ", Px[j]);
		}
		printf("\n!!!!!\n");
#endif // DEBUG

	}

	//decode logic: chose x with the larger P(x|y) to be the decode output bit
	de_interleave(Px_y, Pex);	//use Pex to store the final P(x|y)
	for (int i = 0; i < message_length; ++i) {
#ifdef DEBUG
		printf("%f ", Pex[i]);
#endif // DEBUG

		if (Pex[i] >= 0.5)
			de_message_turbo[i] = 0;
		else
			de_message_turbo[i] = 1;
	}
#ifdef DEBUG
	printf("\n+++++++\n");
#endif // DEBUG

}

/// <summary>
/// main for fun
/// </summary>

//for systematic conv. code test
void main3()
{
	int i;
	float SNR, start, finish;
	long int bit_error, bit_error_soft, bit_error_BCJR, seq, seq_num;
	double BER, BER_soft, BER_BCJR;
	double progress;

	//generate state table
	statetable();

	//random seed
	srand((int)time(0));

	//input the SNR and frame number
	/*printf("\nEnter start SNR: ");
	scanf("%f", &start);
	printf("\nEnter finish SNR: ");
	scanf("%f", &finish);
	printf("\nPlease input the number of message: ");
	scanf("%d", &seq_num);*/
	start = 0;
	finish = 10;
	seq_num = 10000;

	for (SNR = start; SNR <= finish; SNR++)
	{
		//channel noise
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2);

		bit_error = 0;
		bit_error_soft = 0;
		bit_error_BCJR = 0;

		for (seq = 1; seq <= seq_num; seq++)
		{
			//generate binary message randomly
			/****************
			Pay attention that message is appended by 0 whose number is equal to the state of encoder structure.
			****************/


			for (i = 0; i < message_length - state_num; i++)
			{
				message[i] = rand() % 2;
			}
			for (i = message_length - state_num; i < message_length; i++)
			{
				message[i] = 0;		// not need any more for systematic encoder
			}
			//message[0] = 0; message[1] = 0; message[2] = 0; message[3] = 0; message[4] = 0;
			//message[5] = 0; message[6] = 1; message[7] = 1; message[8] = 0; message[9] = 0;


			//convolutional encoder
			encoder_systematic();

			//BPSK modulation
			modulation_systematic();

			//AWGN channel
			channel_systematic();

			//BPSK demodulation, it's needed in hard-decision Viterbi decoder
			demodulation_systematic();

			//convolutional decoder
			BCJR_decoder_systematic(sgm);

			//calculate the number of bit error
			for (i = 0; i < message_length - state_num; i++)
			{
				if (message[i] != de_message_BCJR_systematic[i])
					bit_error++;
			}

			progress = (double)(seq * 100) / (double)seq_num;

			//calculate the BER
			BER = (double)bit_error / (double)((message_length)*seq);

			//print the intermediate result
			printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\r", progress, SNR, bit_error, BER);
			//printf("Progress=%2.1f, SNR=%2.1f, BER=%E, BER_soft=%E, BER_BCJR=%E\r", progress, SNR, BER, BER_soft, BER_BCJR);
		}

		//calculate the BER
		BER = (double)bit_error / (double)((message_length)*seq_num);

		//print the final result
		printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n", progress, SNR, bit_error, BER);
		//printf("SNR = % 2.1f, BER = % E, BER_soft = % E, BER_BCJR = % E\n", SNR, BER, BER_soft, BER_BCJR);
	}
	system("pause");
}

//for turbo code test
void main()
{
	FILE* pFile;
	pFile = fopen("turbo_data.txt", "w");

	int i, SNR_num, SNR_counter;
	float SNR, start, finish, gap;
	long int bit_error, seq, seq_num;
	double BER;
	double progress;

	//generate state table
	statetable();

	//random seed
	srand(0);

	//input the SNR and frame number
	/*printf("\nEnter start SNR: ");
	scanf("%f", &start);
	printf("\nEnter finish SNR: ");
	scanf("%f", &finish);
	printf("\nPlease input the number of message: ");
	scanf("%d", &seq_num);*/
	start = -3;
	finish = 3;
	gap = 0.5;
	seq_num = 1000;
	SNR_num = (finish - start) / gap + 1;
	SNR_counter = 0;
	double* BERs = (double*)malloc(sizeof(double) * SNR_num);

	for (SNR = start; SNR <= finish; SNR+= gap)
	{
		//channel noise
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2);

		bit_error = 0;

		for (seq = 1; seq <= seq_num; seq++)
		{
			//generate binary message randomly
			/****************
			Pay attention that message is appended by 0 whose number is equal to the state of encoder structure.
			****************/


			for (i = 0; i < message_length - state_num; i++)
			{
				message[i] = rand() % 2;
			}
			for (i = message_length - state_num; i < message_length; i++)
			{
				message[i] = 0;		// not need any more for turbo encoder
			}
			//message[0] = 0; message[1] = 0; message[2] = 0; message[3] = 0; message[4] = 0;
			//message[5] = 0; message[6] = 1; message[7] = 1; message[8] = 0; message[9] = 0;


			//convolutional encoder
			encoder_turbo();

			//puncture();

			//BPSK modulation
			//modulation_turbo_puncture();
			modulation_turbo();

			//AWGN channel
			//channel_turbo_puncture();
			channel_turbo();

			//restore_puncture();

			//BPSK demodulation, it's needed in hard-decision Viterbi decoder
			//demodulation_systematic();

			//convolutional decoder
			decoder_turbo(sgm, 18);

			//calculate the number of bit error
			for (i = 0; i < message_length - state_num; i++)
			{
				if (message[i] != de_message_turbo[i])
					bit_error++;
			}

			progress = (double)(seq * 100) / (double)seq_num;

			//calculate the BER
			BER = (double)bit_error / (double)((message_length)*seq);

			//print the intermediate result
			//printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\r", progress, SNR, bit_error, BER);
			//printf("Progress=%2.1f, SNR=%2.1f, BER=%E, BER_soft=%E, BER_BCJR=%E\r", progress, SNR, BER, BER_soft, BER_BCJR);
		}

		//calculate the BER
		BER = (double)bit_error / (double)((message_length)*seq_num);

		//print the final result
		printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n", progress, SNR, bit_error, BER);
		//printf("SNR = % 2.1f, BER = % E, BER_soft = % E, BER_BCJR = % E\n", SNR, BER, BER_soft, BER_BCJR);

		fprintf(pFile, "%f\t", SNR);
		BERs[SNR_counter] = BER;
		SNR_counter++;
	}
	fprintf(pFile, "\n");
	for (int k = 0; k < SNR_num; ++k) {
		fprintf(pFile, "%E\t", BERs[k]);
	}
	fprintf(pFile, "\n");

	fclose(pFile);
	free(BERs);
	system("pause");
}

//for any test use. NOTE: puncture turbo code has unsolved problem!!!!!
int main_test()
{
	for (int i = 0; i < message_length- state_num; i++)
	{
		message[i] = rand()%2;
	}
	for (int i= message_length - state_num; i < message_length; i++)
	{
		message[i] = 0;			// not need in turbo encoder
	}

	/*
	for (int i = 0; i < message_length-2; i++)
	{
		printf("%d ", message[i]);
	}
	printf("----------\n");

	int interleaved_message[message_length-2];
	interleave(message, interleaved_message);
	for (int i = 0; i < message_length-2; i++)
	{
		printf("%d ", interleaved_message[i]);
	}
	printf("----------\n");

	int de_interleaved_message[message_length-2];
	de_interleave(interleaved_message, de_interleaved_message);
	for (int i = 0; i < message_length-2; i++)
	{
		printf("%d ", de_interleaved_message[i]);
	}
	printf("----------\n");
	*/
	encoder_turbo();
#ifdef SHOW_message

	for (int i = 0; i < message_length; ++i) {
		printf("%d ", message[i]);
	}
	printf("\n--------\n");
	for (int i = 0; i < message_length * 3 + state_num; ++i) {
		printf("%d ", codeword_turbo[i]);
	}
	printf("\n--------\n");

	printf("\n***************\n");
#endif // SHOW_message

	N0 = (1.0 / code_rate) / pow(10.0, (float)(16) / 10.0);
	sgm = sqrt(N0 / 2);
	printf("$$$$$$$$$$ %f\n",sgm);
	puncture();
	modulation_turbo_puncture();
	channel_turbo_puncture();
	restore_puncture();
	decoder_turbo(sgm, 1);
	int error_bits = 0;
	for (int i = 0; i < message_length - 2; ++i) {
#ifdef SHOW_message
		printf("%d ", de_message_turbo[i]);

#endif // SHOW_message

		if (de_message_turbo[i] != message[i]) {
			error_bits++;
			printf("error_position:%d\n", i);
		}
	}
#ifdef SHOW_message
	printf("\n----\n");
#endif // SHOW_message

	printf("error bits:%d\n", error_bits);

	return 0;
}
