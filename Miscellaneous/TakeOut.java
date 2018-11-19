//By: Ayush Rautwar
//If you want to substract one from number in brackets in certain situations

import java.util.Scanner;
import java.io.*;


public class TakeOut
{
   public static void main (String [] args) throws Exception
   {
      File file = new File("pleasework.txt");
      Scanner idek= new Scanner(file);
      String boi = idek.nextLine();
      PrintWriter sad = new PrintWriter(new FileOutputStream("hope.txt"));
      while(!boi.equals(null))
      {
         for(int x=0;x<boi.length()-3;x++)
         {
            if(boi.substring(x, x+3).contains("[") && boi.substring(x, x+3).contains("]") && boi.charAt(x+1)>=48 && boi.charAt(x+1)<=57)
               boi = boi.substring(0,x+1)+(Integer.parseInt(boi.substring(x+1,x+2))-1)+boi.substring(x+2);
         }
         for(int x=0;x<boi.length()-5;x++)
         {
         
            if(boi.substring(x, x+5).contains("[") && boi.substring(x, x+5).contains("]") && boi.substring(x, x+5).contains(",") && boi.charAt(x+1)>=48 && boi.charAt(x+1)<=57 && boi.charAt(x+3)>=48 && boi.charAt(x+3)<=57)
               boi = boi.substring(0,x+1)+(Integer.parseInt(boi.substring(x+1,x+2))-1)+","+(Integer.parseInt(boi.substring(x+3,x+4))-1)+boi.substring(x+4);
         }
         sad.print(boi);
         sad.print("\n");
         try{
            boi=idek.nextLine();
         }
         catch(Exception e)
         {
            
            sad.close();
            break;
         }  
      }   
   }
}
   
