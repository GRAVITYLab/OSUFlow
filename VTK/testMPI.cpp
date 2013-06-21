
#include "vtkMPIController.h"

void process(vtkMultiProcessController* controller, void*
vtkNotUsed(arg))
{
        int myId = controller->GetLocalProcessId();
       
        std::cout << "My process id is ";
        std::cout << myId << "." << std::endl;
}

int main( int argc, char* argv[] )
{
        vtkMPIController* controller = vtkMPIController::New();
        controller->Initialize(&argc, &argv);


        controller->SetSingleMethod(process, 0);
        controller->SingleMethodExecute();


        controller->Finalize();
      controller->Delete();
     
      return 0;
}

