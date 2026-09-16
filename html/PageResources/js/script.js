// JavaScript Document
// Menu active Js
// $(document).ready(function () {
//     $(document).on("scroll", onScroll);
    
//     //smoothscroll
//     $('a[href^="#"]').on('click', function (e) {
//         e.preventDefault();
		
// 		$('body').removeClass('has-active-menu');
//     $('#o-wrapper').removeClass('has-slide-left');
//     $('#c-menu--slide-left').removeClass('is-active');
//     $('#c-mask').removeClass('is-active');
	
//         $(document).off("scroll");
        
//         $('a').each(function () {
//             $(this).removeClass('active');
//         })
//         $(this).addClass('active');
		
      
//         var target = this.hash,
//             menu = target;
//         $target = $(target);
//         $('html, body').stop().animate({
//             'scrollTop': $target.offset().top+2
//         }, 1000, 'swing', function () {
//             window.location.hash = target;
//             $(document).on("scroll", onScroll);
//         });
//     });
// });

// function onScroll(event){
//     var scrollPos = $(document).scrollTop();
    
// }

// //STICKY Menu Js
// $(window).scroll(function() {
//     if ($(this).scrollTop() > 1){  
//         $('header').addClass("sticky");
//     }
//     else{
//         $('header').removeClass("sticky");
//     }
// });




//SIde BAr MEnu JS
/**
   * Slide right instantiation and action.
   */
//   var slideRight = new Menu({
//     wrapper: '#o-wrapper',
//     type: 'slide-left',
//     menuOpenerClass: '.c-button',
//     maskId: '#c-mask'
//   });

//   var slideRightBtn = document.querySelector('#c-button--slide-left');
  
//   slideRightBtn.addEventListener('click', function(e) {
//     e.preventDefault;
//     slideRight.open();
//   });
  

	//animation effect(waypoint)
//paste this code under head tag or in a seperate js file.
 // Wait for window load
$(window).load(function() {
  // Animate loader off screen
  $(".se-pre-con").fadeOut("slow");;
  

            function onScrollInit( items, trigger ) {
                items.each( function() {
                var osElement = $(this),
                    osAnimationClass = osElement.attr('data-os-animation'),
                    osAnimationDelay = osElement.attr('data-os-animation-delay');
                  
                    osElement.css({
                        '-webkit-animation-delay':  osAnimationDelay,
                        '-moz-animation-delay':     osAnimationDelay,
                        'animation-delay':          osAnimationDelay
                    });

                    var osTrigger = ( trigger ) ? trigger : osElement;
                    
                    osTrigger.waypoint(function() {
                        osElement.addClass('animated').addClass(osAnimationClass);
                        },{
                            triggerOnce: true,
                            offset: '90%'
                    });
                });
            }

            onScrollInit( $('.os-animation') );
            onScrollInit( $('.staggered-animation'), $('.staggered-animation-container') 
  
  
  
);
});



//Counter JS
 $(document).ready(function() {
jQuery(document).ready(function( $ ) {
        $('.counter').counterUp({
            delay: 10,
            time: 1000
        });
    });
});	

