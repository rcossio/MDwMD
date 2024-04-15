// Page utils
export function clearFormFields() {
    $('form').find('input[type=text], input[type=number], input[type=email], input[type=password], select').val('');
}

// Toast utils
export function showToast(header, message, color = 'yellow') {
    // Check if the toast container exists
    if ($('#submissionToast').length === 0) {
        // Toast HTML structure
        const toastHTML = `
            <div style="position: fixed; top: 80px; right: 20px; min-height: 200px">
                <div style="position: absolute; top: 0; right: 0;">
                    <div id="submissionToast" class="toast" role="alert">
                        <div class="toast-header">
                            <strong class="mr-auto" style="color: black">Submission Status</strong>
                            <button type="button" class="ml-2 mb-1 close" data-dismiss="toast">
                                <span >&times;</span>
                            </button>
                        </div>
                        <div class="toast-body" style="width: 400px;"></div>
                    </div>
                </div>
            </div>
        `;

        // Append the toast HTML to the body
        $('body').append(toastHTML);

        // Initialize the toast (if you're using a library that requires initialization)
        $('#submissionToast').toast({ delay: 3000 }); // Example for Bootstrap toast
    }

    // Set colors based on the 'color' parameter
    const colors = color === 'red' ? 
        {header: 'lightcoral', body: 'lightpink'} : 
        color === 'green' ? 
            {header: 'rgb(173, 234, 173)', body: 'lightcyan'} :
            {header: 'khaki', body: 'khaki'};

    // Update the toast content and colors
    $('#submissionToast .toast-header strong').html(header);
    $('#submissionToast .toast-body').text(message).css('background-color', colors.body);
    $('#submissionToast').find('.toast-header').css('background-color', colors.header);

    // Show the toast
    $('#submissionToast').toast('show');
}


// User utils
function readCookie(name) {
    const nameEQ = name + "=";
    const ca = document.cookie.split(';');
    for(let i=0;i < ca.length;i++) {
        let c = ca[i];
        while (c.charAt(0)==' ') c = c.substring(1,c.length);
        if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
    }
    return null;
}

export function getUserAlias() {
    const userDataCookie = readCookie('userData');
    if (userDataCookie) {
        const userData = JSON.parse(decodeURIComponent(userDataCookie));
        return userData.alias;
    } else {
        return 'Anonymous';
    }
}

export function getUserId() {
    const userDataCookie = readCookie('userData');
    if (userDataCookie) {
        const userData = JSON.parse(decodeURIComponent(userDataCookie));
        return userData._id;
    } else {
        return null;
    }
}

export function clickOnEnter(targetBtnId) {
    document.addEventListener('keydown', function(event) {
        if (event.key === 'Enter') {
            event.preventDefault();
            document.getElementById(targetBtnId).click();
        }
    });
}

export function displayUserAlias(divId) {
    const alias = getUserAlias();
    if (alias === 'Anonymous') {
        document.getElementById(divId).innerHTML = '<a class="btn btn-success" href="/login">Sign In</a>';
    } else {
        document.getElementById(divId).innerHTML = `<i class="fa-solid fa-circle-user"></i> <strong>${alias}</strong>`;
    }
};


export function startSpinnerAnimation(button) {
    button.innerHTML = '<i class="fas fa-spinner fa-spin"></i>';
    button.setAttribute('disabled', true);
}

export function stopSpinnerAnimation(button,text) {
    button.innerHTML = text;
    button.removeAttribute('disabled');
}