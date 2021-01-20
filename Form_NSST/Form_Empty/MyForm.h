#pragma once
#include <Windows.h>
#include <atlstr.h>
#include "Image.h"
#include "Process.h"
#include <msclr\marshal_cppstd.h>
#include <fstream>

#include "ShearParameters.h"
#include "Container.h"
#include "NsstDec.h"
#include "NsstRec.h"

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
	private: System::Windows::Forms::Label^ label2;
	private: System::Windows::Forms::Label^ label3;

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
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
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
			this->menuStrip1->Size = System::Drawing::Size(1050, 28);
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
			this->pictureBox1->Location = System::Drawing::Point(12, 82);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(500, 500);
			this->pictureBox1->SizeMode = System::Windows::Forms::PictureBoxSizeMode::Zoom;
			this->pictureBox1->TabIndex = 1;
			this->pictureBox1->TabStop = false;
			// 
			// pictureBox2
			// 
			this->pictureBox2->Location = System::Drawing::Point(528, 82);
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
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(12, 59);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(99, 17);
			this->label2->TabIndex = 5;
			this->label2->Text = L"Original Image";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(525, 59);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(27, 17);
			this->label3->TabIndex = 6;
			this->label3->Text = L"XR";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1050, 600);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
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
			//long size;
			//BYTE* buffer, *raw_intensity;
			int width, height;
			CString str = openFileDialog1->FileName;

			int displayFlag = 1;

			// LoadBMP can read only 24 bit image depth
			/*buffer = LoadBMP(width, height, size, (LPCTSTR)str);
			raw_intensity = ConvertBMPToIntensity(buffer, width, height);*/

			width = 512;
			height = 512;

			std::ifstream file(str);
			Matrix* image = new Matrix;
			image->CreateMatrix(height, width, 1);

			if (file.is_open()) {
				double temp;
				for (int i = 0; i < height * width; i++) {
					file >> temp;
					image->mat[i] = temp;
				}
			}

			file.close();

			const char* lpfilt = "maxflat";
			ShearParameters shearParameters;

			shearParameters.dcompSize = 4;
			shearParameters.dcomp = new int[4];
			shearParameters.dcomp[0] = 3;
			shearParameters.dcomp[1] = 3;
			shearParameters.dcomp[2] = 4;
			shearParameters.dcomp[3] = 4;

			shearParameters.dsize = new int[4];
			shearParameters.dsize[0] = 32;
			shearParameters.dsize[1] = 32;
			shearParameters.dsize[2] = 16;
			shearParameters.dsize[3] = 16;

			int shearVersion = 0;			// nsst_dec1e
			//int shear_version = 1;		// nsst_dec1
			//int shear_version = 2;		// nsst_dec2


			Cont* dst;
			switch (shearVersion)
			{
			case 0:
				// [dst,shear_f]=nsst_dec1e(image,shear_parameters,lpfilt);
				dst = NsstDec1e(image, shearParameters, lpfilt);			// dst->mats[cellIndex][deepIndex]
				break;

			case 1:
				// [dst, shear_f] = nsst_dec1(image, shear_parameters, lpfilt);
				break;

			case 2:
				// [dst,shear_f]=nsst_dec2(image,shear_parameters,lpfilt);
				break;
			default:
				break;
			}


			Matrix* xr;
			switch (shearVersion)
			{
			case 0:
				xr = NsstRec1(dst, lpfilt);
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



			// display result

			Bitmap^ surface1 = gcnew Bitmap(width, height);
			Bitmap^ surface2 = gcnew Bitmap(width, height);
			pictureBox1->Image = surface1;		// Original image
			pictureBox2->Image = surface2;		// NSST image
			Color c1, c2;

			for (int row = 0; row < height; row++)
				for (int col = 0; col < width; col++) {

					c1 = Color::FromArgb((int)image->mat[row * width + col], (int)image->mat[row * width + col], (int)image->mat[row * width + col]);
					surface1->SetPixel(col, row, c1);

					c2 = Color::FromArgb((int)ceil(xr->mat[row * width + col]), (int)ceil(xr->mat[row * width + col]), (int)ceil(xr->mat[row * width + col]));
					surface2->SetPixel(col, row, c2);
				}


			pictureBox1->Refresh();
			pictureBox2->Refresh();
			//pictureBox1->ImageLocation = openFileDialog1->FileName;
			//labelPath->Text = pictureBox1->ImageLocation;

			/*delete[] buffer;
			delete[] raw_intensity;*/
		}
	}
	};
}
