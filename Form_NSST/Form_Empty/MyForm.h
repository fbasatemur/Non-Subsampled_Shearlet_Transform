#pragma once
#include <Windows.h>
#include <atlstr.h>
#include "Image.h"
#include "Process.h"
#include <msclr\marshal_cppstd.h>

#include "ShearParameters.h"

namespace Form_Empty {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::OpenFileDialog^ openFileDialog1;
	protected:
	private: System::Windows::Forms::MenuStrip^ menuStrip1;
	private: System::Windows::Forms::ToolStripMenuItem^ openToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^ openToolStripMenuItem1;
	private: System::Windows::Forms::PictureBox^ pictureBox1;
	private: System::Windows::Forms::PictureBox^ pictureBox2;
	private: System::Windows::Forms::Label^ label1;
	private: System::Windows::Forms::Label^ labelPath;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container^ components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->openFileDialog1 = (gcnew System::Windows::Forms::OpenFileDialog());
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->openToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->openToolStripMenuItem1 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->pictureBox2 = (gcnew System::Windows::Forms::PictureBox());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->labelPath = (gcnew System::Windows::Forms::Label());
			this->menuStrip1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2))->BeginInit();
			this->SuspendLayout();
			// 
			// openFileDialog1
			// 
			this->openFileDialog1->FileName = L"openFileDialog1";
			// 
			// menuStrip1
			// 
			this->menuStrip1->ImageScalingSize = System::Drawing::Size(20, 20);
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->openToolStripMenuItem });
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Size = System::Drawing::Size(1073, 28);
			this->menuStrip1->TabIndex = 0;
			this->menuStrip1->Text = L"menuStrip1";
			// 
			// openToolStripMenuItem
			// 
			this->openToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->openToolStripMenuItem1 });
			this->openToolStripMenuItem->Name = L"openToolStripMenuItem";
			this->openToolStripMenuItem->Size = System::Drawing::Size(46, 24);
			this->openToolStripMenuItem->Text = L"File";
			// 
			// openToolStripMenuItem1
			// 
			this->openToolStripMenuItem1->Name = L"openToolStripMenuItem1";
			this->openToolStripMenuItem1->Size = System::Drawing::Size(128, 26);
			this->openToolStripMenuItem1->Text = L"Open";
			this->openToolStripMenuItem1->Click += gcnew System::EventHandler(this, &MyForm::openToolStripMenuItem1_Click);
			// 
			// pictureBox1
			// 
			this->pictureBox1->Location = System::Drawing::Point(12, 71);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(500, 500);
			this->pictureBox1->SizeMode = System::Windows::Forms::PictureBoxSizeMode::Zoom;
			this->pictureBox1->TabIndex = 1;
			this->pictureBox1->TabStop = false;
			// 
			// pictureBox2
			// 
			this->pictureBox2->Location = System::Drawing::Point(518, 71);
			this->pictureBox2->Name = L"pictureBox2";
			this->pictureBox2->Size = System::Drawing::Size(500, 500);
			this->pictureBox2->SizeMode = System::Windows::Forms::PictureBoxSizeMode::Zoom;
			this->pictureBox2->TabIndex = 2;
			this->pictureBox2->TabStop = false;
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(12, 32);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(41, 17);
			this->label1->TabIndex = 3;
			this->label1->Text = L"Path:";
			// 
			// labelPath
			// 
			this->labelPath->AutoSize = true;
			this->labelPath->Location = System::Drawing::Point(60, 32);
			this->labelPath->Name = L"labelPath";
			this->labelPath->Size = System::Drawing::Size(20, 17);
			this->labelPath->TabIndex = 4;
			this->labelPath->Text = L"...";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1073, 583);
			this->Controls->Add(this->labelPath);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->pictureBox2);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->menuStrip1);
			this->MainMenuStrip = this->menuStrip1;
			this->Margin = System::Windows::Forms::Padding(4);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox2))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void openToolStripMenuItem1_Click(System::Object^ sender, System::EventArgs^ e) {

		if (openFileDialog1->ShowDialog() == System::Windows::Forms::DialogResult::OK)
		{
			long size;
			int width, height;
			BYTE* buffer, * raw_intensity;
			CString str = openFileDialog1->FileName;


			////////////// this line started develop    

			int displayFlag = 1;

			//BMP image okuma 
			buffer = LoadBMP(width, height, size, (LPCTSTR)str);
			raw_intensity = ConvertBMPToIntensity(buffer, width, height);


			const char* lpfilt = "maxflat";
			ShearParameters shearParameters;
			shearParameters.dcomp = new int[3, 3, 4, 4];
			shearParameters.dsize = new int[32, 32, 16, 16];

			int shearVersion = 0;			// nsst_dec1e
			//int shear_version = 1;		// nsst_dec1
			//int shear_version = 2;		// nsst_dec2

			double** dst, shearF;

			switch (shearVersion)
			{
			case 0:
				//[dst, shear_f] = nsst_dec1e(x, shear_parameters, lpfilt);
				break;

			case 1:
				//[dst, shear_f] = nsst_dec1(x, shear_parameters, lpfilt);
				break;

			case 2:
				//[dst, shear_f] = nsst_dec2(x, shear_parameters, lpfilt);
				break;
			default:
				break;
			}
			
			Bitmap^ surface = gcnew Bitmap(width, height);
			pictureBox2->Image = surface;
			
			/////////////////////  display coefficients - will complete
			//
			// DISPLAY
			// 
			/*if (displayFlag == 1) {

			}*/
			////////////////



			/*
			* FFT IMAGE PROCESS
			*
			*
			// Daha sonra width * height 'luk resmin FFT'si aliniyor
			double* Real_part = new double[width * height];
			double* Im_part = new double[width * height];
			double* image = new double[width * height];

			// width * height 'lik FFT bolgesi secilir
			const int FFTWidth = width;
			const int FFTHeight = height;
			double* cutFFT = new double[FFTWidth * FFTHeight];
			double* fourier255 = new double[FFTWidth * FFTHeight];		// 0-255 fft values
			//double* fourier1 = new double[FFTWidth * FFTHeight];			// 0-1 fft values

			//Gri seviyeyedeki image frekans domeninde ortalanmasi için (-1)^(x+y) ile carpiliyor
			for (int i = 0; i < height; i++)
				for (int j = 0; j < width; j++)
					image[i * width + j] = double(raw_intensity[i * width + j]) * pow(-1, (i + j));


			FFT2D(image, Real_part, Im_part, width, height);

			double deger;
			double maks1 = -1000000000;
			double min1 = 1000000000;

			for (int i = 0; i < height; i++)
				for (int j = 0; j < width; j++) {

					deger = log(0.05 + sqrt(Real_part[i * width + j] * Real_part[i * width + j] + Im_part[i * width + j] * Im_part[i * width + j]));
					cutFFT[i * FFTWidth + j] = deger;

					if (deger > maks1)
						maks1 = deger;
					if (deger < min1)
						min1 = deger;
				}


			for (int i = 0; i < FFTHeight; i++)
				for (int j = 0; j < FFTWidth; j++) {
					fourier255[i * FFTWidth + j] = (cutFFT[i * FFTWidth + j] - min1) / (maks1 - min1) * 255;
					//fourier1[i * FFTWidth + j] = (cutFFT[i * FFTWidth + j] - min1) / (maks1 - min1);
				}

			Color c;
			for (int row = 0; row < FFTHeight; row++)
				for (int col = 0; col < FFTWidth; col++) {

					c = Color::FromArgb(BYTE(fourier255[row * FFTWidth + col]), BYTE(fourier255[row * FFTWidth + col]), BYTE(fourier255[row * FFTWidth + col]));
					surface->SetPixel(col, row, c);
				}

			*/


			double* xr;
			switch (shearVersion)
			{
			case 0:
				//   xr = nsst_rec1(dst, lpfilt);

				break;

			case 1:
				//	xr = nsst_rec1(dst, lpfilt);
				break;

			case 2:
				// xr = nsst_rec2(dst, shear_f, lpfilt);
				break;
			default:
				break;
			}



			pictureBox2->Refresh();
			pictureBox1->ImageLocation = openFileDialog1->FileName;
			labelPath->Text = pictureBox1->ImageLocation;

			delete[] buffer;
			delete[] raw_intensity;
			/*delete[] Real_part;
			delete[] Im_part;
			delete[] fourier255;*/
			//delete[] fourier1;
			//delete[] image;
		}
	}
	};
}
